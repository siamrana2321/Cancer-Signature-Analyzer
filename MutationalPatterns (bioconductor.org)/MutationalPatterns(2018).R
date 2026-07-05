############################################################
# SECTION 1: Load Required Libraries
############################################################

library(MutationalPatterns)
library(BSgenome.Hsapiens.UCSC.hg38)

############################################################
# SECTION 2: Read VCF File
############################################################

vcfs <- read_vcfs_as_granges(
  vcf_files = "D:/Mohammad_Siam_Ahmed_Rana/Mutation_Signature/Resources/sample/SomaticBrca_P1_MP_clean.vcf",
  sample_names = "Tumor_1",
  genome = "hg38",
  predefined_dbs_mbs = TRUE
)

############################################################
# SECTION 3: Generate SBS96 Mutation Matrix
############################################################

mut_mat <- mut_matrix(
  vcf_list = vcfs,
  ref_genome = BSgenome.Hsapiens.UCSC.hg38
)

# Check matrix dimensions
dim(mut_mat)

############################################################
# SECTION 4: Load COSMIC Signatures
############################################################

signatures <- get_known_signatures()

# Check dimensions
dim(signatures)

############################################################
# SECTION 5: Fit COSMIC Signatures
############################################################

fit_res <- fit_to_signatures(
  mut_matrix = mut_mat,
  signatures = signatures
)

############################################################
# SECTION 6: Signature Contributions
############################################################

fit_res$contribution

############################################################
# SECTION 7: Reconstruct Mutation Profile
############################################################

reconstructed <- signatures %*% fit_res$contribution

############################################################
# SECTION 8: Reconstruction Cosine Similarity
############################################################

cosine_similarity <- cos_sim_matrix(
  mut_mat,
  reconstructed
)[1,1]

print(cosine_similarity)

############################################################
# SECTION 9: Correlation
############################################################

correlation <- cor(
  as.vector(mut_mat),
  as.vector(reconstructed)
)

print(correlation)

############################################################
# SECTION 10: RMSE
############################################################

rmse <- sqrt(
  mean(
    (
      as.vector(mut_mat) -
        as.vector(reconstructed)
    )^2
  )
)

print(rmse)

############################################################
# SECTION 11: L1 Norm
############################################################

l1_norm <- sum(
  abs(
    as.vector(mut_mat) -
      as.vector(reconstructed)
  )
)

print(l1_norm)

############################################################
# SECTION 12: L1 Norm Percentage
############################################################

l1_norm_percent <- (
  l1_norm /
    sum(mut_mat)
) * 100

print(l1_norm_percent)

############################################################
# SECTION 13: L2 Norm
############################################################

l2_norm <- sqrt(
  sum(
    (
      as.vector(mut_mat) -
        as.vector(reconstructed)
    )^2
  )
)

print(l2_norm)

############################################################
# SECTION 14: KL Divergence
############################################################

epsilon <- 1e-10

observed_prob <- as.vector(mut_mat) / sum(mut_mat)

reconstructed_prob <- as.vector(reconstructed) /
  sum(reconstructed)

kl_divergence <- sum(
  observed_prob *
    log(
      (observed_prob + epsilon) /
        (reconstructed_prob + epsilon)
    )
)

print(kl_divergence)

############################################################
# SECTION 15: Summary Metrics Table
############################################################

metrics <- data.frame(
  
  Metric = c(
    "Cosine Similarity",
    "Correlation",
    "RMSE",
    "L1 Norm",
    "L1 Norm %",
    "L2 Norm",
    "KL Divergence"
  ),
  
  Value = c(
    cosine_similarity,
    correlation,
    rmse,
    l1_norm,
    l1_norm_percent,
    l2_norm,
    kl_divergence
  )
)

print(metrics)

############################################################
# SECTION 16: Plot SBS96 Profile
############################################################

plot_96_profile(mut_mat)

############################################################
# SECTION 17: Plot Signature Contributions
############################################################

plot_contribution(
  fit_res$contribution
)
