# Export the MetTD model (distance matrix, enzyme table) to binary/CSV for the C kernel.
# One-time job. Output: model_dist.bin (column-major doubles, N*N), enzyme.csv, nodes.csv

suppressWarnings(suppressMessages({
  library(dplyr)
}))

e <- new.env()
load("D:/RSpace/MetTD/human_model.Rda", envir = e)
Model <- get("Model", e)

N <- nrow(Model$distance_martix)
cat("N =", N, "\n")

# enzyme table
enz <- data.frame(idx = Model$enzyme$idx, enzyme = Model$enzyme$Compound.name)
write.csv(enz, "D:/MetTD_recompute/enzyme.csv", row.names = FALSE)

# node table (all 9492 rows, 1-based idx = row/col of matrix)
nodes <- data.frame(idx = 1:N, node = rownames(Model$distance_martix))
write.csv(nodes, "D:/MetTD_recompute/nodes.csv", row.names = FALSE)

cat("dim check:", dim(Model$distance_martix), "\n")
cat("class:", class(Model$distance_martix[[1]]), "\n")

m <- as.matrix(Model$distance_martix)
storage.mode(m) <- "double"
cat("matrix converted, writing binary...\n")
con <- file("D:/MetTD_recompute/model_dist.bin", "wb")
writeBin(as.vector(m), con, size = 8)
close(con)
cat("done. N =", N, " entries =", N * N, "\n")
