# Export ref (node idx + level) for one MetTD dataset, replicating singleProcess.R exactly.
# Usage:
#   Rscript export_ref.R <dir> <out_ref.csv>                       # singleProcess mode (mean diff)
#   Rscript export_ref.R <dir> <out_ref.csv> process_fc <proj> <sx> # process.R mode: rowMeans((m-b)/b),
#                                     model = X{proj}.{1..3}.{sx}, baseline = X{proj}.{1..3}.1
args <- commandArgs(trailingOnly = TRUE)
datadir <- args[1]
outref  <- args[2]
mode    <- if (length(args) >= 3) args[3] else "single"
proj    <- if (length(args) >= 5) as.integer(args[4]) else 0L
sx      <- if (length(args) >= 5) as.integer(args[5]) else 2L

suppressWarnings(suppressMessages({
  library(dplyr)
  library(openxlsx)
}))

local_path <- "D:/RSpace/MetTD/"
source(paste0(local_path, "MetTD.R", sep = ""), local = TRUE)
menv <- new.env(); load(paste0(local_path, "human_model.Rda", sep = ""), envir = menv)
Model <- get("Model", menv)
denv <- new.env(); load(paste0(local_path, "Mixed_db_final.rda", sep = ""), envir = denv)
default_annotation_database <- get("default_annotation_database", denv)

x <- readxl::read_excel(file.path(datadir, "All_annotations.xlsx"))
express <- read.csv(file.path(datadir, "knn_expression_data.csv"))
rownames(express) <- express[[1]]
express <- express[x$variable_id, ]

if (mode == "process_fc") {
  base_cols <- paste0("X", proj, ".", 1:3, ".1")
  model_cols <- paste0("X", proj, ".", 1:3, ".", sx)
  mp <- as.matrix(express[, model_cols])
  bl <- as.matrix(express[, base_cols])
  rel <- (mp - bl) / bl                       # process.R: (model - baseline)/baseline per pair
  lv_fc <- rowMeans(rel)
  lv_md <- rowMeans(mp - bl)
} else {
  # replicate singleProcess.R: model_peak = cols 2..(n+1)/2 (±1.96 SD filter),
  # baseline = remaining cols (±3 SD filter), Level = mp - bl
  modelpk  <- express[, 2:((ncol(express) + 1) / 2)]
  baseline <- express[, (((ncol(express) + 1) / 2) + 1):ncol(express)]
  mp <- vapply(seq_len(nrow(modelpk)), function(i) {
    rdt <- as.numeric(modelpk[i, ]); s <- sd(rdt); m <- mean(rdt)
    mean(rdt[rdt > (m - 1.96 * s) & rdt < (m + 1.96 * s)])
  }, numeric(1))
  bl <- vapply(seq_len(nrow(baseline)), function(i) {
    rdt <- as.numeric(baseline[i, ]); s <- sd(rdt); m <- mean(rdt)
    mean(rdt[rdt > (m - 3 * s) & rdt < (m + 3 * s)])
  }, numeric(1))
  lv_md <- as.numeric(mp - bl)
  lv_fc <- rowMeans((as.matrix(express[, 2:((ncol(express) + 1) / 2)]) -
                     as.matrix(express[, (((ncol(express) + 1) / 2) + 1):ncol(express)])) /
                    as.matrix(express[, (((ncol(express) + 1) / 2) + 1):ncol(express)]))
}
ml <- lv_md
fc <- lv_fc

temp <- x
temp$Level <- ml
ok <- !is.na(ml) & is.finite(ml)
temp <- temp[ok, ]
orig_pos <- which(ok)

ref <- get_ref_idx(temp, Model, default_annotation_database)

# rebuild origin mapping (replicating mutiple_match order) to align fc with ref rows
db_object <- Model$db_object
if (!is.null(db_object$idx)) db_object <- dplyr::select(db_object, -"idx")
colnames(db_object) <- gsub("db_idx", "Lab_ID", colnames(db_object))
meta_db_obj <- match_db(temp, default_annotation_database, by_Synonyms = FALSE)
meta_db_obj <- fill_db_object(temp, meta_db_obj)
meta_db_obj <- dplyr::select(meta_db_obj, -"idx")
colnames(meta_db_obj) <- gsub("db_idx", "Lab_ID", colnames(meta_db_obj))
keep <- !is.na(meta_db_obj$Lab_ID)
kept_i <- which(keep)
mLab <- meta_db_obj$Lab_ID

origin <- c(); ridx <- c()
for (i in seq_along(kept_i)) {
  hits <- which(db_object$Lab_ID == mLab[i])
  if (length(hits) > 0) {
    origin <- c(origin, rep(orig_pos[kept_i[i]], length(hits)))
    ridx <- c(ridx, hits)
  }
}
ref2 <- db_object[ridx, ]
stopifnot(nrow(ref2) == nrow(ref))

fcv <- fc; fcv[!is.finite(fcv)] <- NA
out <- data.frame(
  idx = as.numeric(rownames(ref2)),
  level_meandiff = ref$level,
  level_fc = fcv[origin]
)
write.csv(out, outref, row.names = FALSE)
cat("ref rows:", nrow(out), " unique nodes:", length(unique(out$idx)), "\n")
