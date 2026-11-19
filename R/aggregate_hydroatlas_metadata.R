# =============================================================================
# aggregate_hydroatlas_metadata.R
# -----------------------------------------------------------------------------
# Aggregate HydroATLAS / BasinATLAS v10 hydro-geophysical attributes to a
# per-gage WATERSHED metadata table, sitting alongside the streamflow signatures.
#
# Hybrid aggregation (see docs/plans + CHANGELOG):
#   * Rule 1 (passthrough): attributes that HydroATLAS already accumulates over
#     the full upstream catchment (_u variants) or routes to the pour point
#     (_p, e.g. discharge) are taken AS-IS from the gage's outlet basin (the
#     level-12 subbasin containing the gage). The delineated upstream member set
#     equals the extent HydroATLAS accumulates into that basin's _u fields.
#   * Rule 2 (drop): local _s attributes that have a _u twin are dropped (the
#     _u value from Rule 1 supersedes them).
#   * Rule 3 (manual): _s-only attributes are aggregated across the member
#     basins, SUB_AREA-weighted:
#       - area-weighted mean for continuous / percentage / monthly fields
#         (percentage/areal-extent fields are ALWAYS area-weighted mean,
#          overriding any av/mn/mx suffix; min/max-month temperature too);
#       - spatial min/max ONLY for elevation (a true topographic extreme);
#       - area-weighted majority for categorical class fields
#         (glc/pnv/wet recomputed from the outlet's upstream _u fractions;
#          the rest by area-weighted mode of the _smj code).
#   * Containing-subbasin rule (user decisions 2026-11-19, section 5b): a gage whose
#     containing subbasin is NOT a HydroBASINS headwater and is more than 1.5 x the
#     gage's official (agency) polygon area gets that subbasin ALONE -- its own _s
#     values, NA where HydroATLAS has no subbasin-level value. Record:
#     docs/plans/2026-11-19-hydroatlas-outlet-review/.
#
# This is metadata/ingestion code (R domain). It is NOT a cross-language
# signature and intentionally lives outside rpkg/ (which mirrors Julia
# signatures). No Julia/Python port is required.
#
# Source data (external, D: drive):
#   - basinAt_NorAm_polys.gpkg : HydroBASINS lev12 polygons + 281 BasinATLAS
#     attributes + topology (167,665 NorAm basins, WGS84). HYBAS_ID/NEXT_DOWN
#     are numeric in the file -> coerced to character to match the RDS keys.
#   - upstream_hydrobasins.rds  : list keyed by outlet HYBAS_ID (character);
#     value = member basins (outlet first, then all upstream).
# =============================================================================

suppressMessages({
  library(data.table)
  library(arrow)
})
# sf is needed only by map_gages_to_outlets() (sf:: calls); the containing-subbasin rule runs without it.
if (requireNamespace("sf", quietly = TRUE)) suppressMessages(library(sf))

options(scipen = 999)  # never write HYBAS_IDs in scientific notation

# -----------------------------------------------------------------------------
# HYBAS_ID normalization: numeric -> clean integer-string matching the RDS keys.
# (find_upstream_hydrobasins built the RDS with as.character(); scipen=999 makes
#  as.character() on these <=10-digit doubles emit the full integer string.)
# -----------------------------------------------------------------------------
as_hybas_chr <- function(x) {
  if (is.character(x)) return(x)
  out <- as.character(x)
  out[is.na(x)] <- NA_character_
  out
}

# Canonical gage_id for leading-zero-safe joins: purely-numeric ids have leading
# zeros stripped (so "01011000" == "1011000"); alphanumeric (Canadian) kept.
canonical_gage_id <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^[0-9]+$", x), sub("^0+", "", x), x)
}

# =============================================================================
# 1. ATTRIBUTE CLASSIFIER  ->  rule_table (single source of truth)
# =============================================================================
#' Parse HydroATLAS field names into an aggregation rule table.
#'
#' Naming: <theme3>_<unit2>_<extent><stat>. extent: s=sub-basin, u=upstream,
#' p=pour-point. Returns a data.table with one row per attribute and a `rule`
#' in {passthrough, drop, wmean, min, max, majority}.
classify_hydroatlas_attributes <- function(field_names) {
  topo <- c("HYBAS_ID", "NEXT_DOWN", "NEXT_SINK", "MAIN_BAS", "DIST_SINK",
            "DIST_MAIN", "SUB_AREA", "UP_AREA", "PFAF_ID", "ENDO", "COAST",
            "ORDER_", "SORT", "Shape_Length", "Shape_Area")
  rx <- "^([a-z]{3})_([a-z0-9]{2})_([a-z])([a-z0-9]+)$"
  is_attr <- grepl(rx, field_names) & !(field_names %in% topo)
  attrs <- field_names[is_attr]
  m <- regmatches(attrs, regexec(rx, attrs))
  dt <- data.table(
    attribute = attrs,
    theme  = vapply(m, `[`, "", 2L),
    unit   = vapply(m, `[`, "", 3L),
    extent = vapply(m, `[`, "", 4L),
    stat   = vapply(m, `[`, "", 5L)
  )
  # (theme, unit, stat) family is used to detect whether an _s field has a _u twin
  dt[, family := paste(theme, unit, stat, sep = "_")]
  u_families <- dt[extent == "u", unique(family)]
  dt[, has_u_twin := family %in% u_families]

  dt[, rule := NA_character_]
  dt[extent %in% c("u", "p"), rule := "passthrough"]            # Rule 1
  dt[extent == "s" & has_u_twin, rule := "drop"]                # Rule 2

  # ---- Rule 3 (s-only) ----
  # categorical majority
  dt[extent == "s" & is.na(rule) & unit == "cl", rule := "majority"]
  dt[attribute == "gad_id_smj", rule := "majority"]
  # percentage / areal-extent fields -> ALWAYS area-weighted mean (override mn/mx)
  dt[extent == "s" & is.na(rule) & unit == "pc", rule := "wmean"]
  # monthly (s01..s12) and *av / *yr continuous fields -> area-weighted mean
  dt[extent == "s" & is.na(rule) & grepl("^[0-9]{2}$", stat), rule := "wmean"]
  dt[extent == "s" & is.na(rule) & stat %in% c("av", "yr"), rule := "wmean"]
  # explicit per-field decisions
  dt[attribute %in% c("tmp_dc_smn", "tmp_dc_smx"), rule := "wmean"]  # representative, per decision
  dt[attribute == "ele_mt_smn", rule := "min"]                      # true watershed low point
  dt[attribute == "ele_mt_smx", rule := "max"]                      # true watershed high point

  if (anyNA(dt$rule)) {
    leftover <- dt[is.na(rule), attribute]
    warning("Unclassified HydroATLAS attributes defaulted to area-weighted mean: ",
            paste(leftover, collapse = ", "))
    dt[is.na(rule), rule := "wmean"]
  }
  dt[]
}

# Precompute column groups (and the categorical-from-upstream-fraction map).
.prep_col_groups <- function(rule_table) {
  # glc/pnv have a complete _u fraction basis (all classes), so argmax of the upstream
  # fractions is exact. wet does NOT (wet_pc_u01..u09 only, but source classes reach 11),
  # so wet_cl_smj is handled by the area-weighted mode of its source code instead.
  ufrac <- list(
    glc_cl_smj = sprintf("glc_pc_u%02d", 1:22),
    pnv_cl_smj = sprintf("pnv_pc_u%02d", 1:15)
  )
  maj_all <- rule_table[rule == "majority", attribute]
  pass <- rule_table[rule == "passthrough", attribute]
  # subbasin (_s) twin of each passthrough field, for the containing-subbasin rule; NA where HydroATLAS has
  # no subbasin-level version (lkv/rev _usu sums and the four _p pour-point fields)
  twin <- ifelse(grepl("^[a-z]{3}_[a-z0-9]{2}_u", pass), sub("^([a-z]{3}_[a-z0-9]{2})_u", "\\1_s", pass), NA_character_)
  twin[!twin %in% rule_table$attribute] <- NA_character_
  list(
    pass     = pass,
    pass_twin = stats::setNames(twin, pass),
    wmean    = rule_table[rule == "wmean", attribute],
    min      = rule_table[rule == "min", attribute],
    max      = rule_table[rule == "max", attribute],
    maj_mode = setdiff(maj_all, names(ufrac)),
    maj_ufrac = ufrac[intersect(names(ufrac), maj_all)]
  )
}

# =============================================================================
# 2. GAGE -> OUTLET BASIN  (vectorized point-in-polygon; no file writes)
# =============================================================================
#' Map each gage to the HydroBASINS unit containing it (the watershed outlet).
#' Single st_join of all gage points against the polygon layer (spatial-indexed).
map_gages_to_outlets <- function(metadata_dt, basinAt_NorAm_polys) {
  # HydroBASINS polygons contain minor invalid loops (duplicate vertices) that the
  # s2 spherical engine rejects; GEOS planar point-in-polygon is correct + tolerant.
  old_s2 <- sf::sf_use_s2(FALSE); on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  md <- as.data.table(metadata_dt)[!is.na(latitude) & !is.na(longitude)]
  md[, gage_id := as.character(gage_id)]
  pts <- sf::st_as_sf(md[, .(gage_id, latitude, longitude)],
                      coords = c("longitude", "latitude"), crs = 4326)
  pts <- sf::st_transform(pts, sf::st_crs(basinAt_NorAm_polys))
  polys <- basinAt_NorAm_polys["HYBAS_ID"]
  suppressMessages(
    j <- sf::st_join(pts, polys, join = sf::st_intersects, left = TRUE)
  )
  jdt <- as.data.table(sf::st_drop_geometry(j))
  jdt[, Downstream_HB_ID := as_hybas_chr(HYBAS_ID)]
  jdt <- jdt[, .SD[1L], by = gage_id]          # first match for boundary points
  jdt[is.na(Downstream_HB_ID), Downstream_HB_ID := "NO_BASIN_FOUND"]
  jdt[, .(gage_id, Downstream_HB_ID)]
}

# =============================================================================
# 3. UPSTREAM MEMBER SETS  (cached RDS first; efficient BFS for any missing)
# =============================================================================
#' For each outlet, the set of member basins (outlet + all upstream), using the
#' cached delineation where available and a NEXT_DOWN BFS (adjacency built once)
#' for outlets not yet cached. Equivalent to find_upstream_hydrobasins().
build_upstream_members <- function(outlet_ids, HB_dt, cached = list()) {
  outlet_ids <- unique(outlet_ids)
  # Sanitize the legacy cache: find_upstream_hydrobasins wrote a few member ids in
  # scientific notation (e.g. "7.121e+09"); re-coerce them to clean integer strings so
  # they match the gpkg HYBAS_ID and are not silently dropped at the keyed join.
  if (length(cached)) cached <- lapply(cached, function(v) {
    bad <- !is.na(v) & !grepl("^[0-9]+$", v)
    if (any(bad)) v[bad] <- as_hybas_chr(as.numeric(v[bad]))
    v
  })
  need <- outlet_ids[!outlet_ids %in% names(cached)]
  res <- cached[intersect(outlet_ids, names(cached))]

  if (length(need)) {
    hb <- as_hybas_chr(HB_dt$HYBAS_ID)
    nd <- as_hybas_chr(HB_dt$NEXT_DOWN)
    child_map <- split(hb, nd)                 # downstream id -> upstream children
    max_iter <- length(hb) + 1L
    computed <- vector("list", length(need)); names(computed) <- need
    for (o in need) {
      members <- o; frontier <- o; it <- 0L
      while (length(frontier) && it < max_iter) {
        kids <- unlist(child_map[frontier], use.names = FALSE)
        kids <- kids[!is.na(kids) & !(kids %in% members)]   # tree -> no revisit (defensive)
        if (!length(kids)) break
        members <- c(members, kids); frontier <- kids; it <- it + 1L
      }
      computed[[o]] <- members
    }
    res <- c(res, computed)
  }
  res[outlet_ids[outlet_ids %in% names(res)]]
}

# =============================================================================
# 4. ATTRIBUTE TABLE  (subset to the member union; cache as parquet)
# =============================================================================
load_hydroatlas_attributes <- function(HB_dt, member_ids, attr_cols, cache_path = NULL) {
  dt <- as.data.table(HB_dt)
  dt[, hybas_chr := as_hybas_chr(HYBAS_ID)]
  keep <- intersect(c("hybas_chr", "SUB_AREA", "UP_AREA", "ORDER_", "ENDO", "COAST", attr_cols),
                    names(dt))
  dt <- dt[hybas_chr %in% member_ids, ..keep]
  setkey(dt, hybas_chr)
  if (!is.null(cache_path)) {
    suppressMessages(arrow::write_parquet(dt, cache_path, compression = "snappy"))
  }
  dt[]
}

# =============================================================================
# 5. AGGREGATION  (per UNIQUE outlet; O(1) keyed member slices)
# =============================================================================
#' local_only = TRUE is the containing-subbasin rule (section 5b): the outlet basin ALONE, with every attribute
#' taken from its own subbasin (_s) value; passthrough fields with no _s version are NA.
aggregate_one_outlet <- function(outlet, members, attr_dt, cg, local_only = FALSE) {
  if (local_only) members <- outlet
  sub  <- attr_dt[.(members), nomatch = NULL]
  w    <- as.numeric(sub$SUB_AREA)
  orow <- sub[hybas_chr == outlet]

  res <- list(
    outlet_HYBAS_ID             = outlet,
    n_member_basins             = length(members),
    watershed_area_sum_SUB_AREA = sum(w, na.rm = TRUE),
    outlet_UP_AREA = if (nrow(orow)) as.numeric(orow$UP_AREA) else NA_real_,
    outlet_ORDER   = if (nrow(orow)) as.integer(orow$ORDER_)  else NA_integer_,
    outlet_ENDO    = if (nrow(orow)) as.integer(orow$ENDO)    else NA_integer_,
    outlet_COAST   = if (nrow(orow)) as.integer(orow$COAST)   else NA_integer_
  )

  # Rule 1: passthrough from the outlet basin
  if (length(cg$pass)) {
    if (!nrow(orow)) {
      res[cg$pass] <- NA_real_
    } else if (!local_only) {
      res[cg$pass] <- as.list(orow[, cg$pass, with = FALSE])
    } else {
      # containing-subbasin rule: the basin's own _s twin of each _u field (-9999 NoData -> NA);
      # fields without a subbasin-level version -> NA
      for (pc in cg$pass) {
        tw <- cg$pass_twin[[pc]]
        v <- if (is.na(tw) || !tw %in% names(orow)) NA else orow[[tw]]
        if (!is.na(v) && v == -9999) v <- NA
        res[[pc]] <- v
      }
    }
  }

  # Rule 3a/3b/3c: SUB_AREA-weighted means (NA-aware), vectorized over columns
  if (length(cg$wmean)) {
    M <- as.matrix(sub[, cg$wmean, with = FALSE]) * 1.0
    M[M == -9999] <- NA_real_                           # HydroATLAS NoData sentinel (e.g. sgr_dk_sav)
    if (local_only) {
      wm <- if (nrow(M)) M[1L, ] else rep(NA_real_, length(cg$wmean))   # the basin's own value, exactly
    } else {
      num <- colSums(M * w, na.rm = TRUE)
      den <- colSums((!is.na(M)) * w, na.rm = TRUE)   # na.rm guards NA SUB_AREA weights
      wm  <- num / den
    }
    wm[!is.finite(wm)] <- NA_real_
    res[cg$wmean] <- as.list(wm)
  }

  # Rule 3d/3e: spatial min / max (elevation only)
  for (cm in cg$min) { v <- sub[[cm]]; res[[cm]] <- if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE) }
  for (cm in cg$max) { v <- sub[[cm]]; res[[cm]] <- if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE) }

  # Rule 3f (mode): area-weighted majority of the _smj class code across members.
  # -9999 (HydroATLAS NoData, e.g. wet_cl_smj where no wetland) is excluded so the
  # dominant ACTUAL class wins; NA only if every member is NoData.
  for (cm in cg$maj_mode) {
    dd <- data.table(code = sub[[cm]], a = w)[!is.na(code) & code >= 0 & !is.na(a)]
    res[[cm]] <- if (nrow(dd)) dd[, .(a = sum(a)), by = code][order(-a)]$code[1L] else NA_integer_
  }
  # Rule 3f (upstream fractions): dominant class from the outlet's _u fractions
  # (its own _s fractions under the containing-subbasin rule).
  # Class code is read from the column name (robust if a fraction col is absent).
  for (cm in names(cg$maj_ufrac)) {
    fc0 <- cg$maj_ufrac[[cm]]
    if (local_only) fc0 <- sub("_u([0-9]{2})$", "_s\\1", fc0)
    fcols <- intersect(fc0, names(sub))
    fr <- if (nrow(orow) && length(fcols)) as.numeric(orow[, fcols, with = FALSE]) else NA_real_
    fr[!is.na(fr) & fr == -9999] <- NA_real_
    if (all(is.na(fr)) || max(fr, na.rm = TRUE) <= 0) {   # no class present -> NA (not class 1)
      res[[cm]] <- NA_integer_; res[[paste0(cm, "_frac")]] <- NA_real_
    } else {
      i <- which.max(fr)
      res[[cm]] <- as.integer(sub("^.*_[us]", "", fcols[i]))
      res[[paste0(cm, "_frac")]] <- fr[i]
    }
  }
  res
}

aggregate_hydroatlas_per_outlet <- function(attr_dt, upstream_members, rule_table, outlet_ids) {
  cg <- .prep_col_groups(rule_table)
  if (!identical(data.table::key(attr_dt), "hybas_chr")) setkey(attr_dt, hybas_chr)
  outlet_ids <- unique(outlet_ids)
  rows <- vector("list", length(outlet_ids))
  for (i in seq_along(outlet_ids)) {
    o <- outlet_ids[i]
    m <- upstream_members[[o]]
    if (is.null(m) || !length(m)) next
    rows[[i]] <- aggregate_one_outlet(o, m, attr_dt, cg)
  }
  out <- data.table::rbindlist(rows, fill = TRUE)
  if (nrow(out)) setkey(out, outlet_HYBAS_ID)
  out[]
}

# =============================================================================
# 5b. CONTAINING-SUBBASIN RULE  (user decisions 2026-11-19)
# =============================================================================
#' A gage whose containing level-12 subbasin c is NOT a HydroBASINS headwater (some subbasin has
#' NEXT_DOWN = c) and whose SUB_AREA(c) exceeds `ratio` x the official (agency) polygon area is a small
#' gage inside a through-flow subbasin: it gets c ALONE (aggregate_one_outlet(local_only = TRUE)) instead
#' of c plus everything upstream. Every other gage -- including gages in headwater subbasins -- keeps the
#' upstream aggregation unchanged. Rows are overwritten per GAGE: an outlet can be shared by gages on both
#' sides of the rule. Adds `assignment_method` (upstream / containing_subbasin_only).
#'
#' hb_topology:   data.table(hybas_chr, next_down_chr, SUB_AREA) for every level-12 unit
#' official_area: data.table(canon_id, official_area_km2) for the gages with an agency outline
#'                (geom_area_km2 of the gagesii / wsc_eccc polygons; docs/benchmarks/hydroatlas_assignment_audit.py)
apply_containing_subbasin_rule <- function(gage_dt, attr_dt, rule_table, hb_topology, official_area, ratio = 1.5) {
  g <- data.table::copy(data.table::as.data.table(gage_dt))
  g[, Downstream_HB_ID := as.character(Downstream_HB_ID)]
  cg <- .prep_col_groups(rule_table)
  if (!identical(data.table::key(attr_dt), "hybas_chr")) setkey(attr_dt, hybas_chr)
  stopifnot(all(is.finite(official_area$official_area_km2) & official_area$official_area_km2 > 0),
            !anyDuplicated(canonical_gage_id(official_area$canon_id)))
  oa <- as.numeric(official_area$official_area_km2)[match(canonical_gage_id(g$gage_id),
                                                          canonical_gage_id(official_area$canon_id))]
  ci <- match(g$Downstream_HB_ID, hb_topology$hybas_chr)
  c_sub <- as.numeric(hb_topology$SUB_AREA)[ci]
  headwater <- !(g$Downstream_HB_ID %in% hb_topology$next_down_chr)
  fire <- !is.na(ci) & !headwater & is.finite(oa) & oa > 0 & c_sub > ratio * oa
  idx <- which(fire)
  if (length(idx)) {
    cs <- unique(g$Downstream_HB_ID[idx])
    loc <- data.table::rbindlist(lapply(cs, function(o) aggregate_one_outlet(o, o, attr_dt, cg, local_only = TRUE)),
                                 fill = TRUE)
    m <- match(g$Downstream_HB_ID[idx], loc$outlet_HYBAS_ID)
    stopifnot(!anyNA(m))
    for (col in setdiff(intersect(names(loc), names(g)), "outlet_HYBAS_ID")) {
      v <- loc[[col]][m]
      if (is.integer(g[[col]])) {            # keep the column type (an all-NA column arrives as logical)
        stopifnot(all(is.na(v) | v == round(v)))
        v <- as.integer(v)
      } else if (is.double(g[[col]])) {
        v <- as.numeric(v)
      }
      data.table::set(g, i = idx, j = col, value = v)
    }
    if (all(c("watershed_area_rel_diff", "basin_area", "watershed_area_sum_SUB_AREA") %in% names(g)))
      g[idx, watershed_area_rel_diff := data.table::fifelse(
        is.finite(basin_area) & basin_area > 0,
        abs(watershed_area_sum_SUB_AREA - basin_area) / basin_area, NA_real_)]
  }
  g[, assignment_method := data.table::fifelse(fire, "containing_subbasin_only", "upstream")]
  .order_lead_cols(g)
}

#' Attach `assignment_flag` (ok / topology_defect / no_official_polygon) from the topology audit
#' (docs/benchmarks/hydroatlas_assignment_audit.py), joined on the zero-stripped gage id.
#' The audit must describe THIS table: the same outlet (Downstream_HB_ID) and the same assignment_method for every
#' gage -- a table rebuilt from changed coordinates or metadata needs a fresh audit, never stale flags.
attach_assignment_flags <- function(gage_dt, audit) {
  g <- data.table::copy(data.table::as.data.table(gage_dt))
  if (anyDuplicated(canonical_gage_id(audit$canon_id))) stop("the audit table repeats a gage id")
  i <- match(canonical_gage_id(g$gage_id), canonical_gage_id(audit$canon_id))
  if (anyNA(i)) stop(sum(is.na(i)), " gages are absent from the audit table")
  bad_out <- as.character(audit$Downstream_HB_ID)[i] != as.character(g$Downstream_HB_ID)
  if (any(bad_out)) stop(sum(bad_out), " gages have a different outlet in the audit table: re-run the audit tool on this table")
  if ("assignment_method" %in% names(g)) {
    bad_m <- as.character(audit$assignment_method)[i] != g$assignment_method
    if (any(bad_m)) stop(sum(bad_m), " gages where the R rule and the audit tool disagree")
  }
  f <- as.character(audit$assignment_flag)[i]
  if (anyNA(f)) stop(sum(is.na(f)), " gages have no assignment_flag in the audit table")
  g[, assignment_flag := f]
  .order_lead_cols(g)
}

.lead_cols <- c("gage_id", "latitude", "longitude", "basin_area", "gage_type",
                "Downstream_HB_ID", "n_member_basins", "watershed_area_sum_SUB_AREA",
                "outlet_UP_AREA", "watershed_area_rel_diff", "outlet_ORDER",
                "outlet_ENDO", "outlet_COAST", "assignment_method", "assignment_flag")
.order_lead_cols <- function(dt) {
  lead <- intersect(.lead_cols, names(dt))
  setcolorder(dt, c(lead, setdiff(names(dt), lead)))
  dt[]
}

# =============================================================================
# 6. ASSEMBLE PER-GAGE  +  WRITE  +  DICTIONARY  +  JOIN HELPER
# =============================================================================
assemble_gage_hydroatlas <- function(metadata_dt, gage2outlet, outlet_agg) {
  md <- as.data.table(metadata_dt)
  base_cols <- intersect(c("gage_id", "latitude", "longitude", "basin_area", "gage_type"),
                         names(md))
  md <- md[, ..base_cols]
  md[, gage_id := as.character(gage_id)]
  md <- merge(md, gage2outlet, by = "gage_id", all.x = TRUE)
  out <- merge(md, outlet_agg, by.x = "Downstream_HB_ID", by.y = "outlet_HYBAS_ID",
               all.x = TRUE, sort = FALSE)
  # QA: delineated watershed area vs the gage-reported drainage area (basin_area = truth)
  if (all(c("watershed_area_sum_SUB_AREA", "basin_area") %in% names(out)))
    out[, watershed_area_rel_diff := data.table::fifelse(
      is.finite(basin_area) & basin_area > 0,
      abs(watershed_area_sum_SUB_AREA - basin_area) / basin_area, NA_real_)]
  .order_lead_cols(out)
  out[, gage_id := as.character(gage_id)]
  out[]
}

.hydroatlas_theme_desc <- function() c(
  dis = "Natural discharge (m3/s)", run = "Land-surface runoff (mm)",
  inu = "Inundation extent (%)", lka = "Limnicity / lake area (%)",
  lkv = "Lake volume (1e6 m3)", rev = "Reservoir volume (1e6 m3)",
  dor = "Degree of regulation (%)", ria = "River area (ha)", riv = "River volume (1e3 m3)",
  gwt = "Groundwater table depth (cm)", ele = "Elevation (m)",
  slp = "Terrain slope (deg x10)", sgr = "Stream gradient (dm/km)",
  clz = "Climate zone (class)", cls = "Climate strata (class)",
  tmp = "Air temperature (degC x10)", pre = "Precipitation (mm)",
  pet = "Potential evapotranspiration (mm)", aet = "Actual evapotranspiration (mm)",
  ari = "Global aridity index (x100)", cmi = "Climate moisture index (x100)",
  snw = "Snow cover extent (%)", glc = "Land cover (GLC2000 class / %)",
  pnv = "Potential natural vegetation (class / %)", wet = "Wetland (class / %)",
  `for` = "Forest cover extent (%)", crp = "Cropland extent (%)", pst = "Pasture extent (%)",
  ire = "Irrigated-area extent (%)", gla = "Glacier extent (%)", prm = "Permafrost extent (%)",
  pac = "Protected-area extent (%)", tbi = "Terrestrial biome (class)",
  tec = "Terrestrial ecoregion (class)", fmh = "Freshwater major habitat type (class)",
  fec = "Freshwater ecoregion (class)", cly = "Soil clay fraction (%)",
  slt = "Soil silt fraction (%)", snd = "Soil sand fraction (%)",
  soc = "Soil organic carbon (tonnes/ha)", swc = "Soil water content (%)",
  lit = "Lithological class", kar = "Karst-area extent (%)",
  ero = "Soil erosion (kg/hectare per year; GloSEM v1.2 RUSLE-based estimate for year 2012)", pop = "Population count",
  ppd = "Population density (people/km2)", urb = "Urban extent (%)",
  nli = "Night-time lights (index)", rdd = "Road density (m/km2)",
  hft = "Human footprint index", gad = "GADM administrative area (country id)",
  gdp = "Gross domestic product (USD)",
  hdi = "Human Development Index for 2015 (index value x 1000: range 0-1000, where 1000 = HDI 1.0)"
)

#' Build a data dictionary describing every output column.
build_data_dictionary <- function(rule_table, output_cols) {
  td <- .hydroatlas_theme_desc()
  rule_lab <- c(
    passthrough = "HydroATLAS upstream-accumulated (_u) / pour-point (_p) value at the outlet basin",
    wmean       = "area-weighted mean across upstream basins (SUB_AREA-weighted)",
    min         = "minimum across the aggregated member subbasins",
    max         = "maximum across the aggregated member subbasins",
    majority    = "dominant class: SUB_AREA-weighted mode across upstream member subbasins",
    majority_ufrac = "dominant class: argmax of the outlet basin's upstream (_u) class fractions"
  )
  # what a containing_subbasin_only row holds instead (section 5b)
  local_lab <- "; containing_subbasin_only rows: the containing subbasin's own value"
  diag_desc <- c(
    gage_id = "Gage identifier as published by the agency (USGS site number with its leading zeros, 8 to 15 digits; or Canadian HYDAT station id); join on the zero-stripped form",
    latitude = "Gage latitude (WGS84)", longitude = "Gage longitude (WGS84)",
    basin_area = "Gage-reported drainage area (km2)", gage_type = "USGS or Canada",
    Downstream_HB_ID = "HydroBASINS lev12 outlet basin id containing the gage (read as character; 10-digit)",
    n_member_basins = "Number of HydroBASINS units aggregated (1 for containing_subbasin_only rows)",
    watershed_area_sum_SUB_AREA = "Sum of member SUB_AREA (km2); ~ outlet UP_AREA for upstream rows, the containing subbasin's own area for containing_subbasin_only rows",
    outlet_UP_AREA = "HydroATLAS upstream area at the outlet basin (km2)",
    watershed_area_rel_diff = "|watershed_area_sum_SUB_AREA - gage basin_area| / basin_area (gage-reported drainage area as reference); large values mean the aggregated subbasins differ in size from the gage watershed (see assignment_flag)",
    outlet_ORDER = "Strahler stream order at the outlet basin",
    outlet_ENDO = "Endorheic flag at outlet (0/1/2)", outlet_COAST = "Coastal flag at outlet (0/1)",
    assignment_method = "How the attributes were aggregated: upstream = the containing subbasin plus every subbasin upstream of it; containing_subbasin_only = the containing subbasin alone, the rule applied when that subbasin is not a HydroBASINS headwater and its area exceeds 1.5 x the gage's official polygon area (geom_area_km2 of the agency outline in the unsimplified boundary layer); every attribute is then the subbasin's own value, NA where HydroATLAS has no subbasin-level value",
    assignment_flag = "ok; topology_defect = checked against the agency polygon, the aggregated subbasins either add land outside the watershed (beyond the containing subbasin) amounting to more than half the watershed's area, or miss more than half of the watershed -- except a gage smaller than its containing subbasin whose missed area lies mostly outside the subbasins upstream of it (a small watershed straddling a divide), which is not flagged; no_official_polygon = no agency outline to check against (HydroBASINS-fallback and excluded basins)"
  )
  cg <- .prep_col_groups(rule_table)
  rt <- copy(rule_table)
  dict <- rbindlist(lapply(output_cols, function(col) {
    if (col %in% names(diag_desc)) {
      return(data.table(column = col, aggregation = "diagnostic / key",
                        theme = NA_character_, unit = NA_character_,
                        description = diag_desc[[col]]))
    }
    base_col <- sub("_frac$", "", col)
    r <- rt[attribute == base_col]
    if (!nrow(r)) return(data.table(column = col, aggregation = "derived",
                                    theme = NA_character_, unit = NA_character_,
                                    description = "see HydroATLAS v10 legend"))
    desc <- td[r$theme]; if (is.na(desc)) desc <- "see HydroATLAS v10 legend"
    if (grepl("_frac$", col)) {
      agg <- "areal fraction of the dominant class (0-100)"
    } else if (r$rule == "majority" && col %in% names(cg$maj_ufrac)) {
      agg <- paste0(rule_lab[["majority_ufrac"]],
                    "; containing_subbasin_only rows: argmax of the containing subbasin's own (_s) class fractions")
    } else if (r$rule == "passthrough") {
      tw <- cg$pass_twin[[col]]
      agg <- paste0(rule_lab[["passthrough"]],
                    if (is.na(tw)) "; containing_subbasin_only rows: NA (HydroATLAS has no subbasin-level value)"
                    else sprintf("%s [%s]", local_lab, tw))
    } else {
      agg <- paste0(rule_lab[[r$rule]], local_lab)
    }
    data.table(column = col, aggregation = agg, theme = r$theme, unit = r$unit,
               description = paste0(desc, " [", r$attribute, "]"))
  }), fill = TRUE)
  dict[]
}

#' Write the per-gage metadata (CSV + snappy parquet) and the data dictionary.
write_hydroatlas_metadata <- function(gage_dt, rule_table, out_dir, date_tag) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  csv  <- file.path(out_dir, sprintf("watershed_hydroatlas_metadata_%s.csv", date_tag))
  pq   <- file.path(out_dir, sprintf("watershed_hydroatlas_metadata_%s.parquet", date_tag))
  dict <- file.path(out_dir, "watershed_hydroatlas_metadata_dictionary.csv")
  data.table::fwrite(gage_dt, csv)
  suppressMessages(arrow::write_parquet(gage_dt, pq, compression = "snappy"))
  data.table::fwrite(build_data_dictionary(rule_table, names(gage_dt)), dict)
  message(sprintf("Wrote %d gages x %d cols\n  %s (%.1f MB)\n  %s\n  %s",
                  nrow(gage_dt), ncol(gage_dt), csv, file.size(csv) / 1e6, pq, dict))
  invisible(list(csv = csv, parquet = pq, dictionary = dict))
}

#' Leading-zero-safe join of the HydroATLAS metadata onto a signatures table.
join_hydroatlas_to_signatures <- function(sig_dt, hydroatlas_dt, sig_key = "gage_id") {
  sig <- as.data.table(copy(sig_dt))
  ha  <- as.data.table(copy(hydroatlas_dt))
  sig[, .canon := canonical_gage_id(get(sig_key))]
  ha[,  .canon := canonical_gage_id(gage_id)]
  ha_cols <- setdiff(names(ha), c("gage_id", ".canon"))
  out <- merge(sig, ha[, c(".canon", ha_cols), with = FALSE], by = ".canon",
               all.x = TRUE, sort = FALSE)
  out[, .canon := NULL]
  out[]
}

# =============================================================================
# 7. VALIDATION
# =============================================================================
#' Print and return a validation report for the assembled per-gage metadata.
#' Pass attr_dt + upstream_members to enable the member-coverage and the
#' _u-vs-area-weighted-_s engine cross-check.
validate_hydroatlas_metadata <- function(gage_dt, rule_table,
                                          attr_dt = NULL, upstream_members = NULL,
                                          sample_n = 200L) {
  cat("\n================ HydroATLAS metadata validation ================\n")
  gage_dt <- data.table::as.data.table(data.table::copy(gage_dt))
  if ("Downstream_HB_ID" %in% names(gage_dt)) gage_dt[, Downstream_HB_ID := as.character(Downstream_HB_ID)]
  n <- nrow(gage_dt)
  cat(sprintf("Gages: %d | columns: %d\n", n, ncol(gage_dt)))

  # (5) outlet-resolution failures
  bad <- c("NO_BASIN_FOUND", "NO_HYDROID", "INVALID_HYDROID")
  unresolved <- gage_dt[is.na(Downstream_HB_ID) | Downstream_HB_ID %in% bad |
                          grepl("^ERROR", Downstream_HB_ID)]
  cat(sprintf("Unresolved outlet (no basin / error): %d (%.1f%%)\n",
              nrow(unresolved), 100 * nrow(unresolved) / n))
  resolved <- gage_dt[!gage_id %in% unresolved$gage_id]
  if ("assignment_method" %in% names(gage_dt)) {
    cat(sprintf("assignment_method: %s\n", paste(sprintf("%s=%d", names(table(gage_dt$assignment_method)),
                                                         as.integer(table(gage_dt$assignment_method))), collapse = "  ")))
    if ("assignment_flag" %in% names(gage_dt))
      cat(sprintf("assignment_flag:   %s\n", paste(sprintf("%s=%d", names(table(gage_dt$assignment_flag)),
                                                           as.integer(table(gage_dt$assignment_flag))), collapse = "  ")))
    # checks (1) and (2) compare against the upstream accumulation: containing-subbasin rows are excluded
    resolved <- resolved[assignment_method == "upstream"]
  }

  # (1) area cross-check
  if (all(c("watershed_area_sum_SUB_AREA", "outlet_UP_AREA") %in% names(gage_dt))) {
    ac <- resolved[!is.na(outlet_UP_AREA) & outlet_UP_AREA > 0]
    rel <- abs(ac$watershed_area_sum_SUB_AREA - ac$outlet_UP_AREA) / ac$outlet_UP_AREA
    cat(sprintf("Area check (sum SUB_AREA vs outlet UP_AREA): median %.4f, p95 %.4f, p99 %.4f | >1%%: %d/%d, >5%%: %d\n",
                median(rel, na.rm = TRUE), quantile(rel, .95, na.rm = TRUE),
                quantile(rel, .99, na.rm = TRUE), sum(rel > 0.01, na.rm = TRUE), nrow(ac),
                sum(rel > 0.05, na.rm = TRUE)))
    ord <- order(-rel); k <- min(5L, length(ord))
    cat("  worst area mismatches:\n")
    print(ac[ord[1:k], .(outlet = Downstream_HB_ID, n_member_basins,
                         sum_SUB_AREA = round(watershed_area_sum_SUB_AREA, 1),
                         outlet_UP_AREA = round(outlet_UP_AREA, 1))])
    if ("basin_area" %in% names(ac)) {
      r2 <- ac[!is.na(basin_area) & basin_area > 0, basin_area / outlet_UP_AREA]
      cat(sprintf("  gage basin_area / outlet UP_AREA: median %.3f, IQR [%.3f, %.3f]\n",
                  median(r2, na.rm = TRUE), quantile(r2, .25, na.rm = TRUE),
                  quantile(r2, .75, na.rm = TRUE)))
    }
  }

  # (1b) published QA column: delineated area vs gage-reported basin_area
  if ("watershed_area_rel_diff" %in% names(gage_dt)) {
    rd <- gage_dt[is.finite(watershed_area_rel_diff), watershed_area_rel_diff]
    if (length(rd)) cat(sprintf("watershed_area_rel_diff (sum SUB_AREA vs gage basin_area): median %.3f, p90 %.3f, >25%%: %d/%d\n",
                                median(rd), quantile(rd, .9), sum(rd > 0.25), length(rd)))
  }

  # (3) percentage range -- only fields WE aggregate. Passthrough _u/_p fields
  # mirror HydroATLAS source values (e.g. degree-of-regulation legitimately > 100%).
  pc_cols <- intersect(rule_table[unit == "pc" & rule %in% c("wmean", "min", "max"), attribute],
                       names(gage_dt))
  pc_cols <- c(pc_cols, grep("_frac$", names(gage_dt), value = TRUE))
  if (length(pc_cols)) {
    viol <- vapply(pc_cols, function(cn) {
      v <- gage_dt[[cn]]; sum(v < -0.01 | v > 100.01, na.rm = TRUE)
    }, integer(1))
    cat(sprintf("Percentage fields out of [0,100]: %d offending values across %d cols\n",
                sum(viol), length(pc_cols)))
    if (any(viol > 0)) cat("  offending cols:", paste(names(viol)[viol > 0], collapse = ", "), "\n")
  }

  # (4) per-attribute NA rate (worst offenders)
  attr_out <- intersect(rule_table$attribute, names(gage_dt))
  if (length(attr_out)) {
    nar <- sort(vapply(attr_out, function(cn) mean(is.na(gage_dt[[cn]])), numeric(1)),
                decreasing = TRUE)
    cat(sprintf("Attribute NA rate: median %.3f ; worst 6:\n", median(nar)))
    print(round(head(nar, 6), 3))
  }

  # (6) member basins missing from the attribute table
  if (!is.null(attr_dt) && !is.null(upstream_members)) {
    uni <- unique(unlist(upstream_members, use.names = FALSE))
    miss <- setdiff(uni, attr_dt$hybas_chr)
    cat(sprintf("Member basins missing from attribute table: %d / %d\n", length(miss), length(uni)))
  }

  # (2) _u vs area-weighted _s engine cross-check (elevation)
  if (!is.null(attr_dt) && !is.null(upstream_members) &&
      "ele_mt_uav" %in% names(gage_dt) && "ele_mt_sav" %in% names(attr_dt)) {
    if (!identical(data.table::key(attr_dt), "hybas_chr")) setkey(attr_dt, hybas_chr)
    samp <- head(unique(resolved$Downstream_HB_ID), sample_n)
    man <- vapply(samp, function(o) {
      m <- upstream_members[[o]]; if (is.null(m)) return(NA_real_)
      s <- attr_dt[.(m), nomatch = NULL]; w <- as.numeric(s$SUB_AREA); v <- as.numeric(s$ele_mt_sav)
      sum(v * w, na.rm = TRUE) / sum(w[!is.na(v)], na.rm = TRUE)
    }, numeric(1))
    pass <- resolved[match(samp, Downstream_HB_ID), ele_mt_uav]
    ok <- is.finite(man) & is.finite(pass)
    if (sum(ok) > 5) {
      cat(sprintf("Elevation _u vs area-weighted _s (n=%d): cor=%.4f, median abs diff=%.1f m\n",
                  sum(ok), stats::cor(man[ok], pass[ok]), median(abs(man[ok] - pass[ok]))))
    }
  }
  cat("================================================================\n")
  invisible(list(n = n, n_unresolved = nrow(unresolved)))
}
