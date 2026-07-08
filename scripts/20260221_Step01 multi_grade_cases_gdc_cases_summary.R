suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(data.table)
})

cases_submitter <- c("TCGA-CS-6670","TCGA-DH-5141","TCGA-HT-A619")
endpoint <- "https://api.gdc.cancer.gov/cases"

# grade関連はデータによって入っていたり空だったりするので候補を広めに
fields <- paste(c(
  "case_id","submitter_id","project.project_id",
  "diagnoses.diagnosis_id",
  "diagnoses.primary_diagnosis",
  "diagnoses.morphology",
  "diagnoses.tumor_grade",
  "diagnoses.who_cns_grade",
  "diagnoses.who_nte_grade"
), collapse=",")

query <- list(
  filters = toJSON(list(
    op="in",
    content=list(field="cases.submitter_id", value=cases_submitter)
  ), auto_unbox=TRUE),
  fields = fields,
  expand = "diagnoses",   # ★これが重要
  format = "JSON",
  size = "100"
)

res <- POST(endpoint, body=query, encode="form")
stop_for_status(res)
dat <- fromJSON(content(res, as="text", encoding="UTF-8"), flatten = TRUE)

hits <- as.data.table(dat$data$hits)

cat("Returned columns:\n")
print(names(hits))

grade_cols <- grep("grade|who|cns", names(hits), ignore.case=TRUE, value=TRUE)
cat("\nGrade-related columns found:\n")
print(grade_cols)

keep <- unique(c("submitter_id","project.project_id","diagnoses.primary_diagnosis","diagnoses.morphology", grade_cols))
keep <- intersect(keep, names(hits))
out <- hits[, ..keep]
setorder(out, submitter_id)

cat("\n--- Rows (may include multiple diagnoses per case) ---\n")
print(out)

summ <- out[, lapply(.SD, function(x) paste(unique(na.omit(as.character(x))), collapse=" | ")),
            by=.(submitter_id, project.project_id)]
cat("\n--- Summary ---\n")
print(summ)

fwrite(out,  "multi_grade_cases_cases_expand_diagnoses_raw.csv")
fwrite(summ, "multi_grade_cases_cases_expand_diagnoses_summary.csv")
cat("\nWrote CSVs.\n")