################################################################################
# run_hydroatlas_assignment_rule.R
#
# Apply the containing-subbasin rule (user decisions 2026-11-19) to an EXISTING
# HydroATLAS attribute table -- the deposited one -- without re-running the spatial
# steps (no sf needed). Rows the rule does not change keep every value as stored.
#
#   Rule: a gage whose containing level-12 subbasin is not a HydroBASINS headwater
#   and is more than 1.5 x the gage's official (agency) polygon area gets that
#   subbasin alone; every other gage keeps the upstream aggregation.
#   Record: docs/plans/2026-11-19-hydroatlas-outlet-review/ (PLAN.md v6).
#
# Usage:
#   Rscript run_hydroatlas_assignment_rule.R \
#     --base       <hisss_watershed_hydroatlas_attributes.parquet>   (deposited table) \
#     --basinatlas <basinatlas_attributes.parquet>                   (from the audit tool) \
#     --audit      <hydroatlas_assignment_audit.csv>                 (from the audit tool) \
#     --out        <dir>
# The two audit-tool inputs come from docs/benchmarks/hydroatlas_assignment_audit.py.
# Writes hisss_watershed_hydroatlas_attributes.{csv,parquet} and
# hisss_hydroatlas_attributes_dictionary.csv to --out.
################################################################################

t0 <- Sys.time()
suppressMessages({ library(data.table); library(arrow) })
source("R/aggregate_hydroatlas_metadata.R")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- which(args == flag)
  if (!length(i) || i[1] >= length(args)) stop("missing argument ", flag)
  args[i[1] + 1L]
}
base_path <- get_arg("--base"); ba_path <- get_arg("--basinatlas")
audit_path <- get_arg("--audit"); out_dir <- get_arg("--out")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

base <- as.data.table(arrow::read_parquet(base_path))
base[, gage_id := as.character(gage_id)]
cat(sprintf("base table: %d gages x %d columns (%s)\n", nrow(base), ncol(base), base_path))

ba <- as.data.table(arrow::read_parquet(ba_path))      # HYBAS_ID / NEXT_DOWN stored as strings
stopifnot(is.character(ba$HYBAS_ID), is.character(ba$NEXT_DOWN))
rule_table <- classify_hydroatlas_attributes(names(ba))
hb_topology <- ba[, .(hybas_chr = HYBAS_ID, next_down_chr = NEXT_DOWN, SUB_AREA)]
attr_dt <- load_hydroatlas_attributes(ba, unique(base$Downstream_HB_ID), rule_table$attribute)
rm(ba); invisible(gc())

audit <- fread(audit_path, colClasses = c(canon_id = "character", gage_id = "character",
                                          Downstream_HB_ID = "character"))
official_area <- audit[!is.na(official_area_km2), .(canon_id, official_area_km2)]

gage_dt <- apply_containing_subbasin_rule(base, attr_dt, rule_table, hb_topology, official_area)

# the audit tool decided the rule independently (Python): the two must agree gage for gage
exp_fire <- toupper(as.character(audit$rule_fires))[match(canonical_gage_id(gage_dt$gage_id),
                                                           canonical_gage_id(audit$canon_id))] == "TRUE"
stopifnot(!anyNA(exp_fire))
got_fire <- gage_dt$assignment_method == "containing_subbasin_only"
if (any(exp_fire != got_fire))
  stop(sum(exp_fire != got_fire), " gages where the R rule and the audit tool disagree")
cat(sprintf("rule applied to %d gages (agrees with the audit tool)\n", sum(got_fire)))

gage_dt <- attach_assignment_flags(gage_dt, audit)
stopifnot(nrow(gage_dt) == nrow(base), identical(gage_dt$gage_id, base$gage_id))

validate_hydroatlas_metadata(gage_dt, rule_table)

csv  <- file.path(out_dir, "hisss_watershed_hydroatlas_attributes.csv")
pq   <- file.path(out_dir, "hisss_watershed_hydroatlas_attributes.parquet")
dict <- file.path(out_dir, "hisss_hydroatlas_attributes_dictionary.csv")
data.table::fwrite(gage_dt, csv)
suppressMessages(arrow::write_parquet(gage_dt, pq, compression = "snappy"))
data.table::fwrite(build_data_dictionary(rule_table, names(gage_dt)), dict)
cat(sprintf("wrote %d gages x %d columns\n  %s\n  %s\n  %s\nDONE in %.1f min\n", nrow(gage_dt), ncol(gage_dt),
            csv, pq, dict, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
