library(dada2)

# Arguments
args   <- commandArgs(trailingOnly = TRUE)
indir  <- args[1]
outdir <- args[2]

cat("Input dir:", indir, "\n")
cat("Output dir:", outdir, "\n")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Threads from SLURM (safe fallback): Much better than hardcoding multithread = TRUE — it respects exactly what SLURM allocated and has a safe fallback to 1
n_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Threads:", n_threads, "\n")

# 1. List filtered FASTQ files
fqs <- sort(list.files(indir,
                       pattern = "\\.downsampled\\.fastq\\.gz$",
                       full.names = TRUE))

if(length(fqs) == 0){
    stop("No downsampled FASTQ files found.") # Safety check for empty input: Prevents a silent failure where the script runs but processes nothing
}

sample.names <- sub("\\.downsampled\\.fastq\\.gz$", "", basename(fqs))

cat("Samples found:", length(fqs), "\n")
print(sample.names)

# 2. Learn error model (PacBio-specific)
cat("\n[1/6] Learning error model...\n")
set.seed(42)

err <- learnErrors(fqs,
                   multithread             = n_threads,
                   errorEstimationFunction = PacBioErrfun,
                   BAND_SIZE               = 32,
                   verbose                 = TRUE)

saveRDS(err, file.path(outdir, "error_model.rds"))

pdf(file.path(outdir, "error_plots.pdf"))
plotErrors(err, nominalQ = TRUE)
dev.off()

cat("Error model saved.\n")

# 3. Dereplication
cat("\n[2/6] Dereplication...\n")

derepFs <- derepFastq(fqs, verbose = TRUE) # Explicit dereplication is cleaner, gives you more control, and is the recommended DADA2 workflow
names(derepFs) <- sample.names

# Count input reads
getN <- function(x) sum(getUniques(x))
input_reads <- sapply(derepFs, function(x) sum(x$uniques))

# 4. Denoising
cat("\n[3/6] Denoising...\n")

dadaRes <- dada(derepFs,
                err         = err,
                multithread = n_threads,
                pool        = "pseudo",   # better scaling than TRUE: gives almost identical sensitivity for rare species but scales properly to large sample numbers
                BAND_SIZE   = 32,
                HOMOPOLYMER_GAP_PENALTY = -1) # Specifically important for PacBio — homopolymer regions are the main remaining error source in HiFi reads and this penalises them appropriately during denoising

denoised_reads <- sapply(dadaRes, getN)

# 5. Sequence table + orientation
cat("\n[4/6] Sequence table + orientation...\n")

seqtab <- makeSequenceTable(dadaRes)

# Fix mixed orientation (critical for PacBio without primers)
seqtab <- orientSeqs(seqtab)

saveRDS(seqtab, file.path(outdir, "seqtab_oriented.rds"))

cat("Sequence length distribution:\n")
lens <- nchar(getSequences(seqtab))
print(summary(lens))

pdf(file.path(outdir, "length_distribution.pdf"))
hist(lens, breaks=50, main="ITS length distribution", xlab="Length (bp)") # Very useful QC output — you'll immediately see if your ASVs are the expected ~500–800bp for full ITS, or if something went wrong
dev.off()

# 6. Chimera removal
cat("\n[5/6] Removing chimeras...\n")

seqtab.nochim <- removeBimeraDenovo(seqtab,
                                    method      = "consensus",
                                    multithread = n_threads,
                                    verbose     = TRUE)

saveRDS(seqtab.nochim, file.path(outdir, "seqtab_nochim.rds"))

nonchim_reads <- rowSums(seqtab.nochim)
tabled_reads  <- rowSums(seqtab)

cat("Fraction retained after chimera removal:",
    sum(seqtab.nochim) / sum(seqtab), "\n")

# 7. Read tracking
cat("\n[6/6] Tracking reads...\n")

track <- data.frame(
    sample   = sample.names,
    input    = input_reads,
    denoised = denoised_reads,
    tabled   = tabled_reads,
    nonchim  = nonchim_reads
)

write.csv(track, file.path(outdir, "read_tracking.csv"), row.names = FALSE)

cat("Read tracking summary:\n")
print(track)

# Final summary
cat("\nDADA2 ITS PacBio pipeline complete.\n")
cat("Outputs:\n")
cat(" - error_model.rds\n")
cat(" - seqtab_oriented.rds\n")
cat(" - seqtab_nochim.rds\n")
cat(" - read_tracking.csv\n")
cat(" - error_plots.pdf\n")
cat(" - length_distribution.pdf\n")
