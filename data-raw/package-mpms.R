matU <- rbind(
  c(0.10, 0, 0, 0, 0),
  c(0.05, 0.12, 0.10, 0, 0),
  c(0, 0.35, 0.12, 0.23, 0.12),
  c(0, 0.03, 0.28, 0.52, 0.10),
  c(0, 0, 0.16, 0.11, 0.17)
)

# sexual reproduction component
matF <- rbind(
  c(0, 0, 17.9, 45.6, 0),
  c(0, 0, 0, 0, 0),
  c(0, 0, 0, 0, 0),
  c(0, 0, 0, 0, 0),
  c(0, 0, 0, 0, 0)
)

matA <- matU + matF

popdemo::eigs(matA, what = "lambda")

stages <- c("seed", "small", "medium", "large", "dormant")

colnames(matA) <- stages
rownames(matA) <- stages
colnames(matU) <- stages
rownames(matU) <- stages
colnames(matF) <- stages
rownames(matF) <- stages

mpm1 <- list(matU = matU, matF = matF)


# write to data directory
usethis::use_data(mpm1, overwrite = TRUE)



# # growth/survival component
# matU <- rbind(c(0.1,   0,   0,   0),
#               c(0.6, 0.2, 0.1,   0),
#               c(  0, 0.5, 0.5, 0.1),
#               c(  0,   0, 0.3, 0.8))
#
# # sexual reproduction component
# matF <- rbind(c(  0,   0, 0.2, 0.6),
#               c(  0,   0, 0.1, 0.2),
#               c(  0,   0,   0,   0),
#               c(  0,   0,   0,   0))


# Leslie Matrices
leslie_mpm1 <- list(
  matU = matrix(c(
    0, 0.883, 0, 0, 0, 0, 0, 0, 0, 0, 0.830,
    0, 0, 0, 0, 0, 0, 0, 0, 0.758, 0, 0, 0, 0, 0, 0, 0,
    0, 0.661, 0, 0, 0, 0, 0, 0, 0, 0, 0.540,
    0, 0, 0, 0, 0, 0, 0, 0, 0.398, 0, 0, 0, 0, 0, 0,
    0, 0, 0.253, 0, 0, 0, 0, 0, 0, 0, 0.128
  ), ncol = 8),
  matF = matrix(c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0,
    0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 5,
    0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0,
    0, 0
  ), ncol = 8)
)
usethis::use_data(leslie_mpm1, overwrite = TRUE)
colnames(leslie_mpm1$matU) <- 0:7
leslie_mpm1$matU
colnames(leslie_mpm1$matU) <- 0:7
rownames(leslie_mpm1$matU) <- 0:7
colnames(leslie_mpm1$matF) <- 0:7
rownames(leslie_mpm1$matF) <- 0:7
usethis::use_data(leslie_mpm1, overwrite = TRUE)
