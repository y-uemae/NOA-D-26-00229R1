# GLASS data access & reproduction (Option 2)

The GLASS validation cohort is governed by the Synapse **conditions-for-use** and cannot be redistributed
at the per-sample level. Accordingly, all per-sample GLASS values (expression TPM/log2, ssGSEA Hallmark
scores, ESTIMATE purity, immune scores, and per-sample TP53/IDH variant calls) have been **removed** from
the shipped CSVs, which retain only the GLASS sample identifiers (`case_barcode` / `pair_id`) plus a pointer.

## To reproduce the GLASS figures (Fig. 1B/1C, Fig. 2, Fig. 3, Fig. S7, Fig. S12; Table S5)

1. Create a Synapse account and accept the GLASS conditions-for-use (syn17038081).
2. Download the GLASS WXS/RNA expression and mutation data.
3. Run the cohort-construction pipeline, which regenerates the full per-sample CSVs locally:
   `scripts/…step05c_01 → 02 → 03 → 04`  (→ 05c_glass/glass_final_cohort_*.csv)
   then `step28` (ESTIMATE), `step29` (ssGSEA Hallmark), `step32a` (analysis dataset).
4. Downstream GLASS scripts (step10/11/16/17/32) then run against the regenerated full CSVs.

The redacted identifier-only CSVs let you verify row/sample correspondence after regeneration.
Numeric results (β, CI, p) in the paper were produced from the full CSVs; the audit ledger records the match.
