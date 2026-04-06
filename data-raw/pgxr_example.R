# data-raw/pgxr_example.R
# Generates synthetic paclitaxel dose-response data for pgxR examples
# and vignette. Simulates lung cancer cell lines with varying KRAS, TP53,
# NF1, and WT genomic backgrounds across ancestry groups.
#
# Biological rationale:
#   - NF1-del: sensitized via RAS pathway dysregulation
#   - TP53-mut: intermediate sensitivity, most common NSCLC mutation
#   - WT: intermediate/baseline sensitivity
#   - KRAS-mut: resistant, associated with poor paclitaxel response
#
# Dose range: 0.1-500 nM (clinically relevant paclitaxel concentrations)
# Ancestry: AFR, EUR, EAS, AMR (1KGP/HGDP superpopulation abbreviations)
#
# Run this script to regenerate data/pgxr_example.rda

set.seed(42)

# --- cell line metadata ---
# LL.4 parameters:
#   hill  = b (steepness, negative for increasing response)
#   lower = c (baseline response at zero dose)
#   upper = d (maximum response)
#   ic50  = e (dose at 50% effect, in nM)
cell_lines <- data.frame(
  cell_line = c("pLC1", "pLC2", "pLC3", "pLC4",
                "pLC5", "pLC6", "pLC7", "pLC8"),
  ancestry  = c("AFR", "AFR", "EUR", "EUR",
                "EAS", "EAS", "AMR", "AMR"),
  feature   = c("NF1-del",  "KRAS-mut",
                "NF1-del",  "TP53-mut",
                "TP53-mut", "KRAS-mut",
                "WT",       "WT"),
  hill      = c(-3.5, -1.4, -3.8, -1.8,
                -2.0, -1.5, -2.2, -1.9),
  lower     = c(0, 0, 0, 0, 0, 0, 0, 0),
  upper     = c(98, 72, 98, 87,
                88, 70, 90, 87),
  ic50      = c(2,  120,  3,  22,
                18,  95,  35,  40),
  stringsAsFactors = FALSE
)

# dose points — shifted up so IC50s are well within range
doses <- c(0.5, 1, 2, 5, 10, 25, 100, 500)
n_reps <- 3

# --- LL.4 simulation ---
# f(x) = c + (d - c) / (1 + exp(b * (log(x) - log(e))))
simulate_response <- function(dose, hill, lower, upper, ic50,
                              noise_sd = 4) {
  pred <- lower + (upper - lower) /
    (1 + exp(hill * (log(dose) - log(ic50))))
  resp <- pred + stats::rnorm(length(dose), mean = 0, sd = noise_sd)
  pmax(0, pmin(100, resp))
}

# --- build dataset ---
pgxr_example <- do.call(rbind, lapply(seq_len(nrow(cell_lines)), function(i) {
  cl   <- cell_lines[i, ]
  grid <- expand.grid(
    dose      = doses,
    replicate = seq_len(n_reps)
  )
  grid$response  <- simulate_response(
    grid$dose,
    hill  = cl$hill,
    lower = cl$lower,
    upper = cl$upper,
    ic50  = cl$ic50
  )
  grid$cell_line <- cl$cell_line
  grid$ancestry  <- cl$ancestry
  grid$feature   <- cl$feature
  grid[, c("cell_line", "dose", "response", "ancestry", "feature")]
}))

rownames(pgxr_example) <- NULL

# sanity checks before saving
stopifnot(all(pgxr_example$response >= 0))
stopifnot(all(pgxr_example$response <= 100))
stopifnot(all(c("cell_line", "dose", "response",
                "ancestry", "feature") %in% colnames(pgxr_example)))

cat(sprintf("Generated %d rows across %d cell lines\n",
            nrow(pgxr_example),
            length(unique(pgxr_example$cell_line))))
cat("Cell lines:\n")
print(unique(pgxr_example[, c("cell_line", "ancestry", "feature")]))


usethis::use_data(pgxr_example, overwrite = TRUE)
