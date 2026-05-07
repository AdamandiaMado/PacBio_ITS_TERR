#!/usr/bin/env Rscript

library(dada2)
library(phyloseq)

# ── Threads from SLURM ────────────────────────────────────────────────────────
n_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Threads:", n_threads, "\n")

# ── Paths ─────────────────────────────────────────────────────────────────────
outdir        <- "/data/projects/p605_PerchPilot/PacBio/Roxanne/ITS/raw/filteredRobust/downsampled/dada2"
seqtab_file   <- file.path(outdir, "seqtab_nochim.rds")
unite_db      <- "/data/projects/p605_PerchPilot/PacBio/Colin/ITS/raw/filteredRobust/downsampled/databases/sh_general_release_dynamic_19.02.2025.fasta"
metadata_file <- "/data/projects/p605_PerchPilot/PacBio/Roxanne/ITS/raw/filteredRobust/downsampled/metadata.csv"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ── Checkpoint: skip to export if filtering already done ─────────────────────
if (file.exists(file.path(outdir, "phyloseq_filtered.rds"))) {
    cat("Filtered phyloseq found - skipping to export...\n")
    ps <- readRDS(file.path(outdir, "phyloseq_filtered.rds"))
    cat("Loaded phyloseq object:\n")
    print(ps)

    write.csv(as.data.frame(otu_table(ps)),
              file.path(outdir, "ASV_table_filtered.csv"))
    write.csv(as.data.frame(tax_table(ps)),
              file.path(outdir, "taxonomy_filtered.csv"))
    md <- data.frame(sample_data(ps))
    md$SampleID <- rownames(md)
    write.csv(md,
              file.path(outdir, "metadata_used.csv"),
              row.names = FALSE)

    cat("\n========== SUMMARY ==========\n")
    cat("Final samples:      ", nsamples(ps), "\n")
    cat("Final ASVs:         ", ntaxa(ps), "\n")
    cat("Total reads:        ", sum(sample_sums(ps)), "\n")
    cat("\nSamples per harvest:\n")
    print(table(sample_data(ps)$harvest))
    cat("\nSamples per treatment:\n")
    print(table(sample_data(ps)$treatment))
    cat("\nSamples per sample_type:\n")
    print(table(sample_data(ps)$sample_type))
    cat("\nExport complete.\n")
    quit(save = "no")
}

# ── 1. Load and length-filter sequence table ──────────────────────────────────
cat("\n[1/6] Loading sequence table...\n")
seqtab.nochim <- readRDS(seqtab_file)
cat("Samples:", nrow(seqtab.nochim), "| ASVs:", ncol(seqtab.nochim), "\n")
cat("Total reads:", sum(seqtab.nochim), "\n")

cat("Filtering by sequence length (100-1000bp)...\n")
seq_lengths   <- nchar(colnames(seqtab.nochim))
keep          <- seq_lengths > 100 & seq_lengths < 1000
seqtab.nochim <- seqtab.nochim[, keep]
cat("After length filter:", ncol(seqtab.nochim), "ASVs |",
    sum(seqtab.nochim), "total reads\n")

# ── 2. Assign taxonomy with UNITE ─────────────────────────────────────────────
cat("\n[2/6] Assigning taxonomy (UNITE)...\n")
taxa <- assignTaxonomy(
    seqtab.nochim,
    unite_db,
    multithread = n_threads,
    tryRC       = TRUE,
    verbose     = TRUE
)

saveRDS(taxa, file.path(outdir, "taxa.rds"))
cat("Taxonomy done:", nrow(taxa), "ASVs\n")
cat("Taxonomy columns:", paste(colnames(taxa), collapse = ", "), "\n")
cat("Kingdom distribution:\n")
print(table(taxa[, "Kingdom"], useNA = "ifany"))

write.csv(
    as.data.frame(taxa),
    file.path(outdir, "taxonomy_raw_full.csv")
)

# ── Filter to Fungi BEFORE building phyloseq ──────────────────────────
cat("\nPre-filtering to Fungi only before phyloseq...\n")
fungi_keep <- taxa[, "Kingdom"] %in% c("Fungi", "k__Fungi")
seqtab.nochim <- seqtab.nochim[, fungi_keep]
taxa          <- taxa[fungi_keep, ]
cat("After Fungi pre-filter:", ncol(seqtab.nochim), "ASVs |",
    sum(seqtab.nochim), "total reads\n")

# ── 3. Load and prepare metadata ──────────────────────────────────────────────
cat("\n[3/6] Loading metadata...\n")
metadata <- read.csv(metadata_file, row.names = 1, stringsAsFactors = FALSE)
cat("Metadata rows:", nrow(metadata), "\n")
cat("Columns:", paste(colnames(metadata), collapse = ", "), "\n")

# Keep mesocosm samples and meadow NAT samples; remove blanks and controls
if ("sample_type" %in% colnames(metadata)) {
    metadata <- metadata[metadata$sample_type %in% c("Sample", "Sample_NAT"), ]
    cat("After removing blanks/controls:", nrow(metadata), "samples\n")
}

# Factors
metadata$treatment   <- factor(metadata$treatment)
metadata$harvest     <- factor(metadata$harvest,  levels = c("H1", "H2", "H3"))
metadata$season      <- factor(metadata$season,   levels = c("spring", "summer", "fall"))
metadata$drought     <- factor(metadata$drought)
metadata$temperature <- factor(metadata$temperature)
metadata$heat_waves  <- factor(metadata$heat_waves)
metadata$plate       <- factor(metadata$plate)

# ── 4. Match samples ──────────────────────────────────────────────────────────
cat("\n[4/6] Matching samples...\n")

# Diagnostic to verify name matching before intersect
cat("Example seqtab names:\n")
print(head(rownames(seqtab.nochim)))
cat("Example metadata names:\n")
print(head(rownames(metadata)))

common <- intersect(rownames(seqtab.nochim), rownames(metadata))
cat("Samples in seqtab:  ", nrow(seqtab.nochim), "\n")
cat("Samples in metadata:", nrow(metadata), "\n")
cat("Matching samples:   ", length(common), "\n")

only_seqtab   <- setdiff(rownames(seqtab.nochim), rownames(metadata))
only_metadata <- setdiff(rownames(metadata), rownames(seqtab.nochim))

if (length(only_seqtab) > 0) {
    cat("NOTE - blanks/controls in seqtab but not in metadata (expected):\n")
    print(only_seqtab)
}
if (length(only_metadata) > 0) {
    cat("WARNING - in metadata but NOT in seqtab (dropped during filtering):\n")
    print(only_metadata)
}

if (length(common) == 0) stop("No matching samples! Check sample names.")

seqtab.nochim <- seqtab.nochim[common, ]
metadata       <- metadata[common, , drop = FALSE]

# ── 5. Build phyloseq object ──────────────────────────────────────────────────
cat("\n[5/6] Building phyloseq object...\n")
ps <- phyloseq(
    otu_table(seqtab.nochim, taxa_are_rows = FALSE),
    tax_table(taxa),
    sample_data(metadata)
)
cat("Raw phyloseq object:\n")
print(ps)
cat("Total reads:", sum(sample_sums(ps)), "\n")
saveRDS(ps, file.path(outdir, "phyloseq_raw.rds"))

# ── 6. Filter phyloseq object ─────────────────────────────────────────────────
cat("\n[6/6] Filtering...\n")

# Kingdom filter now redundant (pre-filtered above) but kept as safety net
kingdom_vals <- unique(tax_table(ps)[, "Kingdom"])
cat("Kingdoms detected:", paste(na.omit(kingdom_vals), collapse = ", "), "\n")
if (any(grepl("k__Fungi", kingdom_vals, fixed = TRUE))) {
    ps <- subset_taxa(ps, Kingdom == "k__Fungi")
} else {
    ps <- subset_taxa(ps, Kingdom == "Fungi")
}
cat("After Kingdom == Fungi:", ntaxa(ps), "ASVs |",
    nsamples(ps), "samples |", sum(sample_sums(ps)), "reads\n")

# Remove low abundance ASVs
ps <- prune_taxa(taxa_sums(ps) > 10, ps)
cat("After abundance filter (>10 reads):", ntaxa(ps), "ASVs |",
    sum(sample_sums(ps)), "total reads\n")

# Prevalence filter — orientation-safe using apply()
prev <- apply(otu_table(ps), 2, function(x) sum(x > 0))
ps   <- prune_taxa(prev > 2, ps)
cat("After prevalence filter (>2 samples):", ntaxa(ps), "ASVs |",
    sum(sample_sums(ps)), "total reads\n")

# Remove samples with zero reads after filtering
ps <- prune_samples(sample_sums(ps) > 0, ps)
cat("After removing empty samples:", nsamples(ps), "samples |",
    sum(sample_sums(ps)), "total reads\n")

saveRDS(ps, file.path(outdir, "phyloseq_filtered.rds"))
cat("Filtered phyloseq object:\n")
print(ps)

# ── 7. Export tables ──────────────────────────────────────────────────────────
cat("\nExporting tables...\n")

write.csv(
    as.data.frame(otu_table(ps)),
    file.path(outdir, "ASV_table_filtered.csv")
)
write.csv(
    as.data.frame(tax_table(ps)),
    file.path(outdir, "taxonomy_filtered.csv")
)
md <- data.frame(sample_data(ps))
md$SampleID <- rownames(md)
write.csv(md,
          file.path(outdir, "metadata_used.csv"),
          row.names = FALSE)

# ── 8. Summary ────────────────────────────────────────────────────────────────
cat("\n========== SUMMARY ==========\n")
cat("Final samples:      ", nsamples(ps), "\n")
cat("Final ASVs:         ", ntaxa(ps), "\n")
cat("Total reads:        ", sum(sample_sums(ps)), "\n")
cat("\nSamples per harvest:\n")
print(table(sample_data(ps)$harvest))
cat("\nSamples per treatment:\n")
print(table(sample_data(ps)$treatment))
cat("\nSamples per sample_type:\n")
print(table(sample_data(ps)$sample_type))
cat("\nOutputs saved to:", outdir, "\n")
cat("  phyloseq_raw.rds\n")
cat("  phyloseq_filtered.rds\n")
cat("  taxa.rds\n")
cat("  taxonomy_raw_full.csv\n")
cat("  ASV_table_filtered.csv\n")
cat("  taxonomy_filtered.csv\n")
cat("  metadata_used.csv\n")
cat("\nPipeline complete.\n")
