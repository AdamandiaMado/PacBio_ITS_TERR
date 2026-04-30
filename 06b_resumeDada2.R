library(dada2)

# orientSeqs replacement
orientSeqs_manual <- function(seqtab) {
    rc <- function(sq) {
        chartr("ACGT", "TGCA",
               paste(rev(strsplit(sq, "")[[1]]), collapse=""))
    }

    seqs    <- colnames(seqtab)
    seqs_rc <- sapply(seqs, rc)
    abund   <- colSums(seqtab)

    idx      <- match(seqs_rc, seqs)
    rc_abund <- ifelse(is.na(idx), 0, abund[idx])

    needs_flip <- ifelse(
        !is.na(idx),
        rc_abund >= abund,   # FIX: handle ties
        seqs_rc < seqs
    )

    already_flipped <- rep(FALSE, length(seqs))

    for (i in seq_along(seqs)) {
        if (needs_flip[i] && !already_flipped[i]) {
            colnames(seqtab)[i] <- seqs_rc[i]
            if (!is.na(idx[i])) already_flipped[idx[i]] <- TRUE
        }
    }

    seqtab <- t(rowsum(t(seqtab), colnames(seqtab)))
    return(seqtab)
}
# Paths
indir  <- "/data/projects/p605_PerchPilot/PacBio/Colin/ITS/raw/filteredRobust/downsampled/"
outdir <- "/data/projects/p605_PerchPilot/PacBio/Colin/ITS/raw/filteredRobust/downsampled/dada2"

n_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

# Load saved error model
cat("Loading error model...\n")
err <- readRDS(file.path(outdir, "error_model.rds"))

# List files
fqs <- sort(list.files(indir,
                       pattern    = "\\.fastq\\.gz$",
                       full.names = TRUE))

# Remove any empty files that caused the crash before
fqs <- fqs[file.size(fqs) > 1000]

sample.names <- sub("\\.fastq\\.gz$", "", basename(fqs))
cat("Samples:", length(fqs), "\n")

# Dereplication
cat("\n[1/4] Dereplicating...\n")
derepFs <- derepFastq(fqs, verbose = TRUE)
names(derepFs) <- sample.names

# Denoising
cat("\n[2/4] Denoising...\n")
dadaRes <- dada(derepFs,
                err         = err,
                multithread = n_threads,
                pool        = "pseudo",
                BAND_SIZE   = 32,
                HOMOPOLYMER_GAP_PENALTY = -1)

# Sequence table + orientation
cat("\n[3/4] Sequence table + orientation...\n")
seqtab <- makeSequenceTable(dadaRes)
seqtab <- orientSeqs_manual(seqtab)
saveRDS(seqtab, file.path(outdir, "seqtab_oriented.rds"))

cat("Length distribution:\n")
print(summary(nchar(getSequences(seqtab))))

pdf(file.path(outdir, "length_distribution.pdf"))
hist(nchar(getSequences(seqtab)), breaks = 50,
     main = "ITS length distribution", xlab = "Length (bp)")
dev.off()

# Chimera removal
cat("\n[4/4] Removing chimeras...\n")
seqtab.nochim <- removeBimeraDenovo(seqtab,
                                    method      = "consensus",
                                    multithread = n_threads,
                                    verbose     = TRUE)
saveRDS(seqtab.nochim, file.path(outdir, "seqtab_nochim.rds"))

cat("Fraction retained:", sum(seqtab.nochim) / sum(seqtab), "\n")

# Read tracking
track <- data.frame(
    sample   = sample.names,
    denoised = sapply(dadaRes, function(x) sum(getUniques(x))),
    tabled   = rowSums(seqtab),
    nonchim  = rowSums(seqtab.nochim)
)
write.csv(track, file.path(outdir, "read_tracking.csv"), row.names = FALSE)
cat("Read tracking:\n")
print(track)

cat("\nDone. Resumed from error_model.rds successfully.\n")
