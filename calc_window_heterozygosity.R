
# ----------heterozygosity in windows

genofile <- seqOpen("chr11")

start_pos <- 37000000
end_pos <- 50000000
window_size <- 250000

window_starts <- seq(start_pos, end_pos - window_size, by = window_size)
window_ends <- window_starts + window_size

sample_ids <- seqGetData(genofile, "sample.id")
sample_ids_clean <- gsub("_nodup.bam", "", gsub(".*/", "", sample_ids))

all_positions <- seqGetData(genofile, "position")
all_variant_ids <- seqGetData(genofile, "variant.id")

het_matrix <- matrix(NA, nrow = length(sample_ids_clean), ncol = length(window_starts))
rownames(het_matrix) <- sample_ids_clean
colnames(het_matrix) <- paste0("win_", window_starts)

for (i in seq_along(window_starts)) {
  w_start <- window_starts[i]
  w_end <- window_ends[i]
  win_mask <- (all_positions >= w_start & all_positions < w_end)
  win_var_ids <- all_variant_ids[win_mask]
  
  if (length(win_var_ids) > 0) {
   
    seqSetFilter(genofile, variant.id = win_var_ids, action = "set", verbose = FALSE)
    
    geno <- seqGetData(genofile, "genotype")

    is_het_matrix <- (geno[1, , , drop = FALSE] != geno[2, , , drop = FALSE])
    is_het_matrix <- drop(is_het_matrix)
    is_na_matrix <- (is.na(geno[1, , , drop = FALSE]) | is.na(geno[2, , , drop = FALSE]))
    is_na_matrix <- drop(is_na_matrix)
    
     if (is.vector(is_het_matrix)) {
      is_het_matrix <- matrix(is_het_matrix, nrow = length(sample_ids_clean))
      is_na_matrix <- matrix(is_na_matrix, nrow = length(sample_ids_clean))
    }
  
    is_het_matrix[is_na_matrix] <- NA
    
    het_val <- rowMeans(is_het_matrix, na.rm = TRUE)
    het_matrix[, i] <- het_val
  }
}

data_het_windows <- as.data.frame(het_matrix)
data_het_windows$INDIV <- rownames(data_het_windows)
data_het_windows <- data_het_windows[, c("INDIV", setdiff(names(data_het_windows), "INDIV"))]
