################################################################################
# run_hydroatlas_metadata.R
#
# Build a per-gage WATERSHED hydro-geophysical metadata file from HydroATLAS /
# BasinATLAS v10, sitting alongside the streamflow signatures (keyed by gage_id).
#
# Pipeline: gage point -> outlet HydroBASINS unit (vectorized spatial join) ->
# upstream member set (cached RDS + NEXT_DOWN BFS) -> hybrid per-attribute
# aggregation (prefer HydroATLAS _u; area-weight _s where no _u) -> per gage.
#
# Usage:
#   Rscript run_hydroatlas_metadata.R                 # full run (success gages) -> data_out/
#   Rscript run_hydroatlas_metadata.R --subset 40     # dry-run, 40 spread gages -> test_output/
#   Rscript run_hydroatlas_metadata.R --status all    # all gages with valid coords
#   Rscript run_hydroatlas_metadata.R --metadata <path/to/combined_watershed_metadata.csv>
#   ... --assignment-audit <hydroatlas_assignment_audit.csv>   # REQUIRED for a product build
#   ... --no-assignment-rule                                   # pre-rule table (input to the audit tool)
#
# Containing-subbasin rule (user decisions 2026-11-19; section 5b of R/aggregate_hydroatlas_metadata.R):
# a gage whose containing subbasin is not a HydroBASINS headwater and is more than 1.5 x the gage's
# official polygon area gets that subbasin alone. The official areas and the topology flags come from
# docs/benchmarks/hydroatlas_assignment_audit.py (run it on the --no-assignment-rule output, or on the
# deposited table when the outlets are unchanged). To re-apply the rule to an existing table without the
# spatial steps (no sf), use run_hydroatlas_assignment_rule.R.
#
# Prerequisites: HYDROATLAS_GPKG (config.R), upstream_hydrobasins.rds (repo root),
#   the pinned combined_streamflow_output/combined_watershed_metadata_19feb2026.csv (or --metadata).
################################################################################

t0 <- Sys.time()
suppressMessages({ library(data.table); library(sf); library(arrow) })
source("config.R")
source("R/aggregate_hydroatlas_metadata.R")

# ---- args ----
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) && i[1] < length(args)) args[i[1] + 1L] else default
}
subset_n <- suppressWarnings(as.integer(get_arg("--subset", NA)))
status_filter <- get_arg("--status", "success")
metadata_override <- get_arg("--metadata", NA)
assignment_audit <- get_arg("--assignment-audit", NA)
no_rule <- "--no-assignment-rule" %in% args
if (is.na(assignment_audit) && !no_rule)
  stop("pass --assignment-audit <hydroatlas_assignment_audit.csv> (docs/benchmarks/hydroatlas_assignment_audit.py), ",
       "or --no-assignment-rule for the pre-rule table")

is_dry  <- !is.na(subset_n)
out_dir <- if (is_dry) "test_output" else "data_out"
date_tag <- tolower(format(Sys.Date(), "%d%b%Y"))
if (is_dry) date_tag <- paste0(date_tag, "_dryrun", subset_n)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("========== HYDROATLAS WATERSHED METADATA ==========\n")
cat("Start:", format(t0, "%Y-%m-%d %H:%M:%S"),
    "| mode:", if (is_dry) paste0("DRY-RUN(", subset_n, ")") else "FULL",
    "| status:", status_filter, "| out:", out_dir, "\n")

# ---- locate metadata ----
# PINNED (plan 2026-11-19): the 8,014-gage compile behind every 2026 product; no newest-file selection.
# Pass --metadata explicitly for any other compile (and re-run the audit tool for it).
PINNED_METADATA <- "combined_streamflow_output/combined_watershed_metadata_19feb2026.csv"
metadata_path <- if (!is.na(metadata_override)) metadata_override else PINNED_METADATA
if (!file.exists(metadata_path)) stop("metadata not found: ", metadata_path, " (pass --metadata <path>)")
cat("Metadata:", metadata_path, "\n")
stopifnot(file.exists(metadata_path), file.exists(HYDROATLAS_GPKG))

md <- fread(metadata_path, colClasses = c(gage_id = "character"))
if (status_filter != "all" && "processing_status" %in% names(md)) {
  md <- md[processing_status == status_filter]
}
md <- md[!is.na(latitude) & !is.na(longitude)]
if (is_dry && nrow(md) > subset_n) md <- md[round(seq(1, .N, length.out = subset_n))]  # spread for variety
cat("Gages to process:", nrow(md), "\n")

# ---- load HydroATLAS gpkg once (peak memory: ~1.5-2 GB) ----
cat("Loading gpkg (167k polygons; this is the slow step)...\n")
basinAt <- sf::st_read(HYDROATLAS_GPKG, quiet = TRUE)

cat("Mapping gages to outlet basins (single st_join)...\n")
gage2outlet <- map_gages_to_outlets(md, basinAt)

# strip geometry -> attribute data.table, then free the heavy sf object
HB_strip <- basinAt; sf::st_geometry(HB_strip) <- NULL
HB_dt <- data.table::as.data.table(HB_strip)
rm(HB_strip, basinAt); gc()

rule_table <- classify_hydroatlas_attributes(names(HB_dt))
cat("Attribute rules:",
    paste(sprintf("%s=%d", names(table(rule_table$rule)), as.integer(table(rule_table$rule))),
          collapse = "  "), "\n")

outlet_ids <- gage2outlet[grepl("^[0-9]+$", Downstream_HB_ID), unique(Downstream_HB_ID)]
cat(sprintf("Resolved outlets: %d unique (of %d gages; %d unresolved)\n",
            length(outlet_ids), nrow(gage2outlet),
            gage2outlet[Downstream_HB_ID == "NO_BASIN_FOUND", .N]))

# ---- upstream member sets (cache first, BFS the rest) ----
cat("Building upstream member sets...\n")
cached <- if (file.exists("upstream_hydrobasins.rds")) readRDS("upstream_hydrobasins.rds") else list()
upstream_members <- build_upstream_members(outlet_ids, HB_dt, cached)
member_union <- unique(unlist(upstream_members, use.names = FALSE))
cat(sprintf("  outlets: %d | unique member basins: %d | from cache: %d\n",
            length(upstream_members), length(member_union), sum(outlet_ids %in% names(cached))))

# ---- attribute table (subset to members; cache parquet) ----
attr_cache <- file.path(out_dir, sprintf("hydroatlas_member_attrs_%s.parquet", date_tag))
attr_dt <- load_hydroatlas_attributes(HB_dt, member_union, rule_table$attribute, cache_path = attr_cache)
hb_topology <- HB_dt[, .(hybas_chr = as_hybas_chr(HYBAS_ID), next_down_chr = as_hybas_chr(NEXT_DOWN), SUB_AREA)]
rm(HB_dt); gc()
cat(sprintf("Attribute table: %d basins x %d cols  (cached: %s)\n",
            nrow(attr_dt), ncol(attr_dt), attr_cache))

# ---- aggregate (per unique outlet) + assemble (per gage) ----
cat("Aggregating per outlet...\n")
outlet_agg <- aggregate_hydroatlas_per_outlet(attr_dt, upstream_members, rule_table, outlet_ids)
gage_dt <- assemble_gage_hydroatlas(md, gage2outlet, outlet_agg)
if (!no_rule) {
  audit <- fread(assignment_audit, colClasses = c(canon_id = "character", gage_id = "character",
                                                  Downstream_HB_ID = "character"))
  gage_dt <- apply_containing_subbasin_rule(gage_dt, attr_dt, rule_table, hb_topology,
                                            audit[!is.na(official_area_km2), .(canon_id, official_area_km2)])
  gage_dt <- attach_assignment_flags(gage_dt, audit)
}

# ---- validate + write ----
validate_hydroatlas_metadata(gage_dt, rule_table, attr_dt, upstream_members)
write_hydroatlas_metadata(gage_dt, rule_table, out_dir, date_tag)

cat(sprintf("\nDONE in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
