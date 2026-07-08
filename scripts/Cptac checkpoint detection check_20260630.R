###############################################################################
## CPTAC-GBM : protein-level detection of checkpoint genes  (CORRECTED)
## Identifier column = 'Composite.Element.REF' with 'SYMBOL|SYMBOL' format.
## Confirms abundant proteins are quantified while LAG3 / checkpoint genes
## are below the MS detection limit -> evidence for the text response to
## Reviewer #1 C1 / Reviewer #2 C2.
###############################################################################
suppressPackageStartupMessages(library(data.table))

results_root <- here::here("results", "TP53", "20260630")
out_dir      <- file.path(results_root, "CPTAC_LAG3_protein_output")
prot_file    <- file.path(results_root,
                          "gbm_cptac_2021_extracted/gbm_cptac_2021/data_protein_quantification.txt")

prot <- fread(prot_file, check.names = FALSE, data.table = FALSE)

# Gene symbol = part before '|' in Composite.Element.REF
id_col <- "Composite.Element.REF"
symbols <- sub("\\|.*$", "", as.character(prot[[id_col]]))
sample_cols <- setdiff(names(prot), id_col)

detect <- function(sym) {
  i <- which(symbols == sym)
  if (length(i) == 0) return(c(in_matrix = 0, n_detected = 0, detection_rate = 0))
  v <- suppressWarnings(as.numeric(unlist(prot[i[1], sample_cols])))
  c(in_matrix = 1, n_detected = sum(!is.na(v)),
    detection_rate = round(mean(!is.na(v)), 3))
}

checkpoint <- c("LAG3","PDCD1","CD274","PDCD1LG2","CTLA4","TIGIT","HAVCR2")
controls   <- c("EGFR","GFAP","VIM","B2M","GAPDH","ACTB")

tab <- data.frame(
  gene  = c(checkpoint, controls),
  panel = c(rep("checkpoint", length(checkpoint)),
            rep("positive_control", length(controls))),
  t(sapply(c(checkpoint, controls), detect)),
  row.names = NULL
)

cat("Total proteins quantified in global proteome :", nrow(prot), "\n")
cat("Samples                                       :", length(sample_cols), "\n\n")
print(tab, row.names = FALSE)

n_cp <- sum(tab$panel == "checkpoint" & tab$in_matrix == 1)
n_pc <- sum(tab$panel == "positive_control" & tab$in_matrix == 1)
cat(sprintf("\nCheckpoint genes present in proteome : %d of 7\n", n_cp))
cat(sprintf("Positive controls present            : %d of %d\n",
            n_pc, length(controls)))
cat(sprintf("LAG3 present in proteome             : %s\n",
            ifelse("LAG3" %in% symbols, "YES", "NO")))

fwrite(tab, file.path(out_dir, "checkpoint_protein_detection_CORRECTED.csv"))
cat("\nWrote checkpoint_protein_detection_CORRECTED.csv to:\n", out_dir, "\n")
