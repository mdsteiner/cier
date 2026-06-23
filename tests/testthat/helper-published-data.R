# Shared helpers for the published-results reproduction suite
# (test-published-results-table.R, test-published-results-reproduction.R).
#
# The slow tier re-fetches the external datasets (fetch-only; never bundled) and
# recomputes every cier cell of inst/extdata/published-results.csv with the same
# conventions and seeds the vignette uses. A deliberate white-box check: it pins
# CSV<->code agreement; faithfulness to the papers is pinned by the frozen
# transcription and the allowed-differs list in the fast tier.

pubres_seed <- 1L            # RPR / RIR / Gnormed resampling seed
pubres_draw_seed <- 20260610L  # Schroeders 1,000-test-draw loop seed

pubres_urls <- list(
  battery_zip = "https://polybox.ethz.ch/index.php/s/88QNfcknKt9JHQm/download",
  rpr_zip     = "https://polybox.ethz.ch/index.php/s/OUeUYS9ZHkyu7YB/download",
  schroeders_resp = "https://osf.io/download/brfmz/"
)

pubres_cache <- function() {
  dir <- file.path(tools::R_user_dir("cier", "cache"), "published-results")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}

# Download into the cache unless already present; NULL (never an error) on
# failure so callers can skip cleanly when offline.
pubres_fetch <- function(url, file) {
  dest <- file.path(pubres_cache(), file)
  if (file.exists(dest) && file.size(dest) > 0) return(dest)
  old <- options(timeout = 600)
  on.exit(options(old), add = TRUE)
  ok <- tryCatch(utils::download.file(url, dest, mode = "wb", quiet = TRUE),
                 error = function(e) -1L, warning = function(w) -1L)
  if (identical(ok, 0L)) return(dest)
  if (file.exists(dest)) unlink(dest)
  NULL
}

pubres_unzip <- function(zip, subdir, probe) {
  out <- file.path(pubres_cache(), subdir)
  if (!file.exists(file.path(out, probe))) {
    utils::unzip(zip, exdir = out)
  }
  file.path(out, probe)
}

# Skip (never fail) when an upstream file is missing or has silently changed: a
# re-uploaded dataset must be re-inspected, not silently compared.
pubres_pin <- function(actual, expected, what) {
  if (!identical(actual, expected)) {
    testthat::skip(sprintf("upstream data changed: %s is %s, pinned %s",
                           what, paste(actual, collapse = "/"),
                           paste(expected, collapse = "/")))
  }
  invisible(TRUE)
}

# ---- the committed table ----------------------------------------------------

pubres_table <- function() {
  path <- system.file("extdata", "published-results.csv", package = "cier")
  if (identical(path, "")) {
    stop("published-results.csv is missing from the package", call. = FALSE)
  }
  classes <- c(paper = "character", study = "character",
               design = "character", level = "character",
               index = "character", paper_statistic = "character",
               paper_value = "numeric", cier_value = "numeric",
               delta = "numeric", status = "character", n = "integer",
               note = "character")
  utils::read.csv(path, stringsAsFactors = FALSE, colClasses = classes)
}

# Digest of the paper-side (transcription) columns only: the cier-side columns are
# free to change on a re-knit, the transcription is not.
pubres_paper_digest <- function(tab) {
  key <- sprintf("%s|%s|%s|%s|%s|%s|%s|%s", tab$paper, tab$study, tab$design,
                 tab$level, tab$index, tab$paper_statistic,
                 ifelse(is.na(tab$paper_value), "NA",
                        formatC(tab$paper_value, format = "f", digits = 4)),
                 ifelse(is.na(tab$n), "NA", as.character(tab$n)))
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  # radix sort = C-locale byte order; writeBin = no platform EOL translation
  payload <- paste0(paste(sort(key, method = "radix"), collapse = "\n"), "\n")
  writeBin(charToRaw(payload), tmp)
  unname(tools::md5sum(tmp))
}

pubres_cell <- function(tab, ...) {
  sel <- rep(TRUE, nrow(tab))
  for (cond in list(...)) sel <- sel & tab[[cond[1]]] == cond[2]
  tab[sel, , drop = FALSE]
}

# ---- shared statistics ------------------------------------------------------

# Rank (Mann-Whitney) AUC of `score` for label 1 over label 0; the score must already
# be oriented careless-high. Non-finite scores and NA labels drop.
pubres_auc <- function(score, label) {
  ok <- is.finite(score) & !is.na(label)
  s <- score[ok]
  l <- as.logical(label[ok])
  pos <- s[l]
  neg <- s[!l]
  if (!length(pos) || !length(neg)) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}

# Empirical sensitivity at 95% specificity: cutoff = 95th percentile of the
# careful group's (careless-high) score; sensitivity = share of careless above.
pubres_sen95 <- function(score, label) {
  ok <- is.finite(score) & !is.na(label)
  s <- score[ok]
  l <- as.logical(label[ok])
  if (!sum(l) || !sum(!l)) return(NA_real_)
  mean(s[l] > stats::quantile(s[!l], 0.95, type = 7))
}

pubres_orient <- function(v, dir) if (dir == "lower") -as.numeric(v) else as.numeric(v)

# Score every respondent on a fixed pair list (Meade & Craig per-person correlation
# across stacked pairs) -- used for the careful-reference synonym / antonym companion
# cells, where pair discovery and scoring run on different row sets so the wrapper's
# whole-sample pairing does not apply.
pubres_score_pairs <- function(m, pairs) {
  if (is.null(pairs) || !nrow(pairs)) return(rep(NA_real_, nrow(m)))
  a <- m[, pairs$item_i, drop = FALSE]
  b <- m[, pairs$item_j, drop = FALSE]
  vapply(seq_len(nrow(m)), function(i) {
    x <- a[i, ]
    y <- b[i, ]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3 || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) {
      return(NA_real_)
    }
    stats::cor(x[ok], y[ok])
  }, numeric(1))
}

# ---- Goldammer battery (2024) ----------------------------------------------

pubres_battery_dir <- function() {
  zip <- pubres_fetch(pubres_urls$battery_zip, "goldammer-battery.zip")
  if (is.null(zip)) testthat::skip("Goldammer battery data not fetchable")
  pubres_unzip(zip, "battery", "Careless responding revisted")
}

pubres_battery_read <- function(dir, file) {
  utils::read.csv(file.path(dir, file), check.names = FALSE,
                  stringsAsFactors = FALSE)
}

# Raw matrix + facet-level item metadata + per-respondent total seconds for one wave
# of BFI-2 items. `suffix` is "" (Studies 1/2) or "_t1" / "_t2" (Studies 4/5 wide
# files). The total is the row sum of the per-item Qualtrics page-submit timers
# shipped in the same files (untimed cells drop; a respondent with no timed item gets
# NA and abstains downstream).
pubres_battery_wave <- function(d, suffix = "") {
  base <- sub(paste0(suffix, "$"), "",
              grep(sprintf("^[eacno]_[a-z]+[0-9]+%s$", suffix), names(d),
                   value = TRUE))
  m <- as.matrix(d[, paste0(base, suffix)])
  storage.mode(m) <- "double"
  colnames(m) <- base
  twin <- vapply(base, function(x) {
    paste0(x, "R", suffix) %in% names(d) || paste0(x, suffix, "R") %in% names(d)
  }, logical(1))
  tm <- d[, paste0(base, suffix, "_Page_Submit")]
  tm[] <- lapply(tm, function(x) as.numeric(gsub(",", ".", as.character(x))))
  tm <- as.matrix(tm)
  tm[!is.finite(tm) | tm <= 0] <- NA_real_
  seconds <- rowSums(tm, na.rm = TRUE)
  seconds[rowSums(!is.na(tm)) == 0L] <- NA_real_
  list(m = m, seconds = seconds,
       items = data.frame(scale = sub("[0-9]+$", "", base),
                          reverse_keyed = unname(twin),
                          max = 6L))
}

# All per-respondent index values for one wave: the five matrix-only indices,
# facet-split RPR, per-domain pooled Gnormed / Ht, whole-sample synonym / antonym
# scores at the paper's critical_r = .60, and the total completion time (fast totals
# careless, so direction lower).
pubres_battery_values <- function(w) {
  pooled <- function(fn) {
    dom <- substr(w$items$scale, 1, 1)
    out <- matrix(NA_real_, nrow(w$m), length(unique(dom)))
    for (k in seq_along(unique(dom))) {
      sel <- dom == unique(dom)[k]
      im <- w$items[sel, , drop = FALSE]
      im$scale <- substr(im$scale, 1, 1)
      # only typed cier errors (e.g. a backend span limit) abstain a domain; anything
      # else is a defect and propagates
      na_domain <- function(e) rep(NA_real_, nrow(w$m))
      out[, k] <- tryCatch(suppressWarnings(fn(w$m[, sel, drop = FALSE], im)),
                           cier_error = na_domain)
    }
    g <- rowMeans(out, na.rm = TRUE)
    g[is.nan(g)] <- NA_real_
    g
  }
  # no tryCatch on the pairing wrappers: a zero-pair sample must ABSTAIN (NA values),
  # never error -- an error here is a real defect and must fail the test loudly
  rpr <- cier_personal_reliability(w$m, w$items, seed = pubres_seed)$value
  vals <- list()
  vals$RPR <- list(v = rpr, dir = "upper")
  vals$MD <- list(v = suppressWarnings(cier_mahalanobis(w$m)$value), dir = "upper")
  vals$r_pbis <- list(v = cier_person_total(w$m)$value, dir = "lower")
  vals$Longstring <- list(v = cier_longstring(w$m)$value, dir = "upper")
  vals$IRV <- list(v = cier_irv(w$m)$value, dir = "lower")
  vals$Time <- list(v = cier_total_time(w$seconds)$value, dir = "lower")
  vals$Gnormed <- list(v = pooled(function(m, im) {
    cier_gnormed(m, im, seed = pubres_seed)$value
  }), dir = "upper")
  vals$Ht <- list(v = pooled(function(m, im) cier_ht(m, im)$value),
                  dir = "lower")
  syn <- suppressWarnings(cier_psychsyn(w$m, critical_r = 0.60)$value)
  ant <- suppressWarnings(cier_psychant(w$m, critical_r = 0.60)$value)
  vals$Synonyms <- list(v = syn, dir = "lower")
  vals$Antonyms <- list(v = ant, dir = "upper")
  vals
}

# Careful-reference companion scores (Studies 4/5): pairs discovered on the careful
# conditions' second wave, scored on every respondent of both waves.
pubres_battery_companions <- function(m1, m2, careful_rows) {
  ref <- m2[careful_rows, , drop = FALSE]
  syn_pairs <- cier_synonym_pairs(ref, critical_r = 0.60)
  ant_pairs <- cier_synonym_pairs(ref, critical_r = 0.60, antonym = TRUE)
  comp <- list()
  comp[["Synonyms (careful-reference)"]] <-
    list(t1 = pubres_score_pairs(m1, syn_pairs),
         t2 = pubres_score_pairs(m2, syn_pairs), dir = "lower")
  comp[["Antonyms (careful-reference)"]] <-
    list(t1 = pubres_score_pairs(m1, ant_pairs),
         t2 = pubres_score_pairs(m2, ant_pairs), dir = "upper")
  comp
}

# All recomputed battery cells as rows (study, design, level, index, value).
pubres_battery_row <- function(study, design, level, index, value) {
  data.frame(study = study, design = design, level = level, index = index,
             cier_value = round(value, 3), stringsAsFactors = FALSE)
}

pubres_battery_cells_s12 <- function(dir, study) {
  d <- pubres_battery_read(dir, sprintf("Study %s/Study_%s.csv", study, study))
  vals <- pubres_battery_values(pubres_battery_wave(d))
  rows <- list()
  for (lv in c("partial", "full")) {
    labcol <- if (lv == "partial") "careless_33" else "careless_100"
    lab <- as.integer(d[[labcol]])
    for (nm in names(vals)) {
      score <- pubres_orient(vals[[nm]]$v, vals[[nm]]$dir)
      rows[[length(rows) + 1L]] <-
        pubres_battery_row(study, "BS", lv, nm, pubres_auc(score, lab))
    }
  }
  do.call(rbind, rows)
}

pubres_battery_cells_s45 <- function(dir, study) {
  d <- pubres_battery_read(dir,
                           sprintf("Study %s/Study_%s_wide.csv", study, study))
  w1 <- pubres_battery_wave(d, "_t1")
  w2 <- pubres_battery_wave(d, "_t2")
  v1 <- pubres_battery_values(w1)
  v2 <- pubres_battery_values(w2)
  cond <- as.integer(d$condition)
  comp <- pubres_battery_companions(w1$m, w2$m, which(cond %in% c(1L, 2L)))
  score_of <- function(nm, wave) {
    if (nm %in% names(v1)) {
      list(v = (if (wave == 1) v1 else v2)[[nm]]$v, dir = v1[[nm]]$dir)
    } else {
      list(v = comp[[nm]][[paste0("t", wave)]], dir = comp[[nm]]$dir)
    }
  }
  rows <- list()
  for (lv in c("partial", "full")) {
    cc <- if (lv == "full") 3L else 4L
    lab_bs <- rep(NA_integer_, nrow(d))
    lab_bs[cond %in% c(1L, 2L)] <- 0L
    lab_bs[cond == cc] <- 1L
    rws <- which(cond == cc)
    lab_ws <- c(rep(0L, length(rws)), rep(1L, length(rws)))
    for (nm in c(names(v1), names(comp))) {
      s1 <- score_of(nm, 1)
      s2 <- score_of(nm, 2)
      auc_bs <- pubres_auc(pubres_orient(s2$v, s2$dir), lab_bs)
      auc_ws <- pubres_auc(pubres_orient(c(s1$v[rws], s2$v[rws]), s1$dir),
                           lab_ws)
      rows[[length(rows) + 1L]] <- pubres_battery_row(study, "BS", lv, nm,
                                                      auc_bs)
      rows[[length(rows) + 1L]] <- pubres_battery_row(study, "WS", lv, nm,
                                                      auc_ws)
    }
  }
  do.call(rbind, rows)
}

pubres_battery_cells <- function(dir) {
  rbind(pubres_battery_cells_s12(dir, "1"), pubres_battery_cells_s12(dir, "2"),
        pubres_battery_cells_s45(dir, "4"), pubres_battery_cells_s45(dir, "5"))
}

# ---- Goldammer RPR paper (2024) ----------------------------------------------

pubres_rpr_facets <- function() {
  list(mlq_is = paste0("mlq_is", 1:4),   mlq_iib = paste0("mlq_iib", 1:4),
       mlq_iia = paste0("mlq_iia", 1:4), mlq_ic = paste0("mlq_ic", 1:4),
       mlq_im = paste0("mlq_im", 1:4),   mlq_app = paste0("mlq_app", 1:4),
       mlq_mbp = paste0("mlq_mbp", 1:4), mlq_lf = paste0("mlq_lf", 1:4),
       alq_t = paste0("alq_t", 1:5),     alq_m = paste0("alq_m", 1:4),
       alq_pb = paste0("alq_pb", 1:3),   alq_sa = paste0("alq_sa", 1:4),
       lmx_a = paste0("lmx_a", 1:3),     lmx_l = paste0("lmx_l", 1:3),
       lmx_pr = paste0("lmx_pr", 1:3),   lmx_e = paste0("lmx_e", 1:3),
       oca = paste0("oca", 1:5),         occ = paste0("occ", 1:4),
       ocn = paste0("ocn", 1:5),
       ocb_h = paste0("ocb_h", 1:5),     ocb_g = paste0("ocb_g", 1:5),
       ocb_u = paste0("ocb_u", 1:5),     ocb_e = paste0("ocb_e", 1:5))
}

pubres_rpr_data <- function() {
  zip <- pubres_fetch(pubres_urls$rpr_zip, "goldammer-rpr.zip")
  if (is.null(zip)) testthat::skip("Goldammer RPR data not fetchable")
  base <- "Comparing conventional and resampled personal reliability"
  csv <- pubres_unzip(zip, "rpr",
                      file.path(base, "real_data_example_240118.csv"))
  utils::read.csv(csv, check.names = FALSE, stringsAsFactors = FALSE)
}

# Recomputed Table 9 cells: PR = even-odd consistency (the paper's definition), RPR
# at 25 / 50 / 100 resamples, each on the 6- / 16- / 23-facet pools.
pubres_rpr_cells <- function(d) {
  facets <- pubres_rpr_facets()
  pools <- list(`6 facets` = names(facets)[1:6],
                `16 facets` = names(facets)[1:16],
                `23 facets` = names(facets))
  rev_items <- c("oca2", paste0("ocb_u", 1:5))
  lab <- as.integer(d$carel_grp)
  rows <- list()
  for (pool in names(pools)) {
    fs <- pools[[pool]]
    its <- unlist(facets[fs], use.names = FALSE)
    m <- as.matrix(d[, its])
    storage.mode(m) <- "double"
    im <- data.frame(scale = rep(fs, lengths(facets[fs])),
                     reverse_keyed = its %in% rev_items, max = 5L)
    vals <- list(PR = cier_even_odd(m, im)$value)
    for (nr in c(25L, 50L, 100L)) {
      vals[[paste0("RPR", nr)]] <-
        cier_personal_reliability(m, im, n_resamples = nr,
                                  seed = pubres_seed)$value
    }
    for (nm in names(vals)) {
      rows[[length(rows) + 1L]] <-
        data.frame(level = pool, index = nm,
                   AUC = round(pubres_auc(vals[[nm]], lab), 3),
                   SEN95 = round(pubres_sen95(vals[[nm]], lab), 3),
                   stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

# ---- Schroeders et al. (2022) -----------------------------------------------

# The authors' scale-blocked item order for even-odd consistency (their analysis
# script 06_online_trad: items_hh / _em / _ex / _ag / _co / _op).
pubres_schroeders_eo_items <- function() {
  blk <- function(ids) paste0("HE01_", ids)
  c(blk(c("06", "30_r", "54", "12_r", "36", "60_r", "18", "42_r", "24_r", "48_r")),
    blk(c("05", "29", "53_r", "11", "35_r", "17", "41_r", "23", "47", "59_r")),
    blk(c("04", "28_r", "52_r", "10_r", "34", "58", "16", "40", "22", "46_r")),
    blk(c("03", "27", "09_r", "33", "51", "15_r", "39", "57_r", "21_r", "45")),
    blk(c("02", "26_r", "08", "32_r", "14_r", "38", "50", "20_r", "44_r", "56_r")),
    blk(c("01_r", "25", "07", "31_r", "13", "37", "49_r", "19_r", "43", "55_r")))
}

# The authors' binarisation rules (their OSF scripts 03 + 06), one tiny function per
# rule so the fast tier can pin each boundary deterministically. IRV: the paper
# computed both tails and kept the better-performing side; on this data that is the
# LOWER tail (instructed-inattentive respondents straightline), which reproduces all
# five Table-3 metrics -- the upper tail (the only variant left uncommented in their
# script) does not.
pubres_flag_longstring <- function(v) v >= 6
pubres_flag_irv <- function(v) v <= stats::quantile(v, 0.10, type = 7)
pubres_flag_antonyms <- function(v) {
  v[is.na(v)] <- 0
  v > 0
}
pubres_flag_evenodd <- function(v) v > 0 & !is.na(v)

pubres_schroeders_data <- function() {
  csv <- pubres_fetch(pubres_urls$schroeders_resp, "schroeders-data_mod_resp.csv")
  if (is.null(csv)) testthat::skip("Schroeders data not fetchable")
  d <- utils::read.csv2(csv, stringsAsFactors = FALSE)
  d$X <- NULL
  d
}

# Replicate the authors' evaluation loop: 1,000 random test samples of 162 regular +
# 18 careless respondents; per draw, each index is computed within the draw and
# binarised with the authors' published rule; the reported cell is the mean metric
# across draws.
pubres_schroeders_cells <- function(d, iterations = 1000L) {
  he <- grep("^HE01_", names(d), value = TRUE)
  recoded <- as.matrix(d[, he])
  storage.mode(recoded) <- "double"
  unrec <- recoded
  rcols <- grep("_r$", colnames(unrec))
  unrec[, rcols] <- 6 - unrec[, rcols]
  qorder <- order(as.integer(sub("^HE01_(\\d+).*$", "\\1", colnames(unrec))))
  unrec <- unrec[, qorder]
  eo_items <- pubres_schroeders_eo_items()
  eo_meta <- data.frame(scale = rep(c("hh", "em", "ex", "ag", "co", "op"),
                                    each = 10),
                        reverse_keyed = FALSE, max = 5L)
  careless <- which(d$Careless == 1L)
  regular <- which(d$Careless == 0L)
  flags <- function(rows) {
    m <- recoded[rows, , drop = FALSE]
    eo <- cier_even_odd(m[, eo_items, drop = FALSE], eo_meta)$value
    ls <- cier_longstring(unrec[rows, , drop = FALSE])$value
    fl <- list()
    fl$Mahalanobis <- cier_mahalanobis(m, alpha = 0.05)$flagged %in% TRUE
    fl$Antonyms <-
      pubres_flag_antonyms(cier_psychant(m, critical_r = 0.20)$value)
    fl[["Even-odd"]] <- pubres_flag_evenodd(eo)
    fl$Longstring <- pubres_flag_longstring(ls)
    fl$IRV <- pubres_flag_irv(cier_irv(m)$value)
    fl
  }
  metric_names <- c("accuracy", "sensitivity", "specificity", "precision",
                    "balanced accuracy")
  acc <- array(NA_real_, c(iterations, 5L, 5L),
               dimnames = list(NULL, names(flags(seq_len(20))), metric_names))
  set.seed(pubres_draw_seed)
  for (i in seq_len(iterations)) {
    rows <- c(sample(regular, 162L), sample(careless, 18L))
    l <- d$Careless[rows] == 1L
    fl <- flags(rows)
    for (nm in names(fl)) {
      f <- fl[[nm]]
      sens <- sum(f & l) / sum(l)
      spec <- sum(!f & !l) / sum(!l)
      acc[i, nm, ] <- c(mean(f == l), sens, spec,
                        if (sum(f)) sum(f & l) / sum(f) else NA_real_,
                        (sens + spec) / 2)
    }
  }
  out <- expand.grid(index = dimnames(acc)[[2]], paper_statistic = metric_names,
                     stringsAsFactors = FALSE)
  out$cier_value <- round(mapply(function(ix, st) {
    mean(acc[, ix, st], na.rm = TRUE)
  }, out$index, out$paper_statistic), 3)
  out
}

# ---- Bruhlmann et al. (2020) -------------------------------------------------

# Flagged counts at the authors' literal cutoffs, on the bundled BFI-44 plus the two
# bundled attention checks. Their pass/fail rules treat an undefined (NA) consistency
# value as flagged; the check pass-sets follow their analysis script (bogus fails at
# agreement >= 3, instructed fails off the directed 0). A missing check abstains
# (cier_attention returns NA) and is not a failure, matching Bruhlmann's
# missing-is-pass rule -- hence na.rm on the check counts.
pubres_bruhlmann_cells <- function() {
  bfi <- as.matrix(cier::bfi_careless[, 1:44])
  items <- data.frame(scale = sub("^v_BFI_([A-Za-z]+)[0-9].*$", "\\1",
                                  colnames(bfi)),
                      reverse_keyed = grepl("_R$", colnames(bfi)),
                      max = 5L)
  ls <- cier_longstring(bfi)$value
  oec <- cier_even_odd(bfi, items)$value
  rir <- cier_personal_reliability(bfi, items, n_resamples = 100L,
                                   seed = pubres_seed)$value
  ptc <- cier_person_total(bfi)$value
  bog <- cier_attention(cier::bfi_careless["v_Bogus_Item"],
                        pass = list(c(1, 2)))$value
  iri <- cier_attention(cier::bfi_careless["v_IRI"], pass = list(0))$value
  data.frame(index = c("Longstring", "Odd-even consistency",
                       "Resampled individual reliability",
                       "Person-total correlation", "Bogus item",
                       "Instructed response item"),
             cier_value = c(sum(ls > 22),
                            sum(oec > 0 | is.na(oec)),
                            sum(rir > 0 | is.na(rir)),
                            sum(ptc < 0 | is.na(ptc)),
                            sum(bog >= 1, na.rm = TRUE),
                            sum(iri >= 1, na.rm = TRUE)),
             stringsAsFactors = FALSE)
}
