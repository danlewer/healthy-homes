# This code estimates power for a matched cohort study comparing homes that are renovated with homes that are not yet renovated.
# This is exposure density sampling at the household level, with outcome aggregated by household
# It assumes:
#  - 3 types of homes for exact matching (eg. EPC rating)
#  - 10% risk of moving house each year
#  - controls are sampled without replacement, but re-released each year to the next batch of retrofits

# ---------
# libraries
# ---------

library(data.table)
library(MASS)
library(PASSED)

# ------
# inputs
# ------

# distribution of home sizes
# https://www.gov.uk/government/statistical-data-sets/social-and-private-renters
probs <- data.frame(hh_size = 1:6, prob = c(0.452, 0.224, 0.136, 0.098, 0.057, 0.033))
sum(probs$hh_size * probs$prob) # mean = 2.18

# groupings for exact matching

home_types <- data.frame(home_type = 1:3, prob = c(1, 1, 1))

# number of retrofits per year

rf <- data.table(
  year = 2022:2029,
  n = c(100, 130, 300, 1000, 2500, 2536, 2536, 1898)
)

# example dispersion

nb_disp <- function (mu, size) c('1+' = 1 - pnbinom(0, size = size, mu = mu), '95' = qnbinom(0.95, size = size, mu = mu), '99' = qnbinom(0.99, size = size, mu = mu))
nb_disp(mu = 1, size = 0.5)

# -------------
# simulate data
# -------------

sim <- function (n_homes = 11000, 
                 prob_new_tenancy = 0.1, 
                 study_end = Inf,
                 vuln_prev = 1,
                 outcome_rate = 0.3,
                 outcome_disp = 0.5,
                 effect_irr = 0.9,
                 printdim = F) {
  
  n_years <- length(rf$year)
  types <- sample(home_types$home_type, n_homes, T, home_types$prob)
  dl <- data.table(
    home_id = rep(seq_len(n_homes), each = n_years),
    type = rep(types, each = n_years),
    year = rep(rf$year, n_homes)
  )
  dl$new_tenancy <- rbinom(nrow(dl), 1, prob_new_tenancy)
  dl$tenancy_id <- cumsum((dl$year == rf$year[1]) | (dl$new_tenancy == 1))
  
  # family size & vulnerability
  hh <- data.table(
    tenancy_id = seq_len(max(dl$tenancy_id)),
    hh_size = sample(probs$hh_size, max(dl$tenancy_id), T, probs$prob)
  )
  hh$vuln_n <- rbinom(nrow(hh), hh$hh_size, vuln_prev)
  dl <- dl[hh, on = 'tenancy_id ']
  
  # retrofit sequence
  eligible_homes <- seq_len(n_homes)
  retrofit_sequence <- NULL
  for (i in rf$year) {
    e <- data.table(retrofit_year = i, home_id = sample(eligible_homes, rf[year == i, n], F))
    retrofit_sequence <- rbind(retrofit_sequence, e)
    eligible_homes <- eligible_homes[!(eligible_homes %in% retrofit_sequence$home_id)]
  }
  dl <- dl[retrofit_sequence, on = 'home_id']
  dl$retrofit <- dl$year >= dl$retrofit_year
  
  # limit to those with vulnerability
  dl <- dl[vuln_n > 0]
  
  # add end of follow-up
  dl <- dl[year <= study_end]
  
  # do matching - NOTE CONTROLS ARE SAMPLED WITHOUT REPLACEMENT
  f <- function (i = 2025, j = 1, d = dl) {
    i1 <- d[retrofit_year == i & year == i & type == j, tenancy_id]
    c1 <- d[year == i & retrofit == F & type == j]
    sampled_controls <- min(length(i1), nrow(c1))
    controls <- c1[sample(seq_along(home_id), sampled_controls, F), tenancy_id]
    i1 <- sample(i1, sampled_controls, F)
    x <- rbind(d[tenancy_id %in% i1 & year >= i], d[tenancy_id %in% controls & retrofit == F & year >= i])
    x <- x[, .(years = .N), c('home_id', 'type', 'tenancy_id', 'retrofit')]
    hh <- unique(d[, c('tenancy_id', 'vuln_n')])
    hh[x, on = 'tenancy_id']
  }
  inputs <- CJ(types = unique(dl$type), years = rf$year)
  m <- mapply(f, i = inputs$years, j = inputs$types, SIMPLIFY = F)
  m <- rbindlist(m)
  m$py <- m$vuln_n * m$years

  # model outcomes
  if (printdim) print (dim(m))
  m$target <- ifelse(m$retrofit, outcome_rate * effect_irr, outcome_rate) * m$py
  m$outcome <- rnbinom(nrow(m), size = 1, mu = m$target)
  
  # fit NB model
  mod <- glm.nb(outcome ~ retrofit + offset(log(py)), data = m)
  summary(mod)$coef[2,c(1, 4)]
}

# --------------
# power function
# --------------

power <- function(B = 10, inputs, effects = seq(0.7, 1, 0.05), alpha = 0.05) {
  n <- length(effects)
  i <- data.frame(lapply(inputs, rep, times = n))
  i$effect_irr <- rep(effects, each = nrow(inputs))
  r <- lapply(1:B, function (x) {
    print (paste0(x, '/', B))
    cbind(i, t(do.call(mapply, c(FUN = sim, i[, -1]))))
  })
  r <- rbindlist(r)
  r$B <- rep(1:B, each = nrow(i))
  ns <- names(r)[!(names(r) %in% c('Estimate', 'Pr(>|z|)', 'B'))]
  r[, .(B = .N, power = mean(`Pr(>|z|)` < alpha), b = exp(mean(Estimate))), by = ns]
}

# --------------------
# varying outcome rate
# --------------------

inputs <- expand.grid(scenario = NA,
                      outcome_rate = c(0.05, c(1:5/10)),
                      vuln_prev = c(1, 0.1),
                      study_end = 2028)
inputs$outcome_disp <- ifelse(inputs$vuln_prev == 1, 0.5, 5)
inputs$scenario = LETTERS[1:nrow(inputs)]

# this takes a long time to run, so results are pre-calculated and loaded
# results <- power(B = 100, inputs = inputs, effects = seq(0.7, 1, 0.025))
results <- readRDS(url("https://github.com/danlewer/healthy-homes/raw/main/power3_results.RDS"))

pf <- function (scenarios, cols, d) {
  plot(1, type = 'n', xlim = c(0.7, 1), ylim = c(0, 1), axes = F, xlab = 'Effect size (IRR)', ylab = 'Power')
  axis(1, seq(0.7, 1, 0.1), pos = 0)
  axis(2, 0:5/5, pos = 0.7, las = 2)
  rect(0.7, 0, 1, 1)
  segments(0.7, 0.05, x1 = 1, lty = 3)
  lapply(seq_along(scenarios), function (i) {
    a <- d[scenario == scenarios[i]]
    m <- glm(power ~ poly(effect_irr, 2), data = a, family = 'binomial')
    nd <- data.frame(effect_irr = seq(0.7, 1, 0.001))
    nd$p <- predict(m, newdata = nd, type = 'response')
    p90 <- nd[which.min(abs(nd$p - 0.9)),'effect_irr']
    print(p90)
    points(a$effect_irr, a$power, col = cols[i])
    lines(nd$effect_irr, nd$p, col = cols[i])
    segments(0.7, 0.9, x1 = p90, lty = 3, col = cols[i])
    segments(p90, 0, y1 = 0.9, lty = 3, col = cols[i])
  })
}

png('hh_power_17aug2024.png', height = 6, width = 12, units = 'in', res = 300)
cols <- brewer.pal(6, 'Set2')
par(mfrow = c(1, 2), mar = c(5, 5, 3, 0), oma = c(0, 0, 0, 10), xpd = NA)
pf(LETTERS[1:6], cols = cols, d = results)
title(main = 'Primary analysis', line = 1)
pf(LETTERS[7:12], col = cols, d = results)
title(main = 'Vulnerable subgroup', line = 1)
ys <- seq(0.5, 0.9, length.out = 6)
segments(1.05, ys, 1.08, col = cols)
text(1.103, ys, c(0.05, c(1:5/10)), adj = 0)
text(1.05, 1, 'Mean events ppy', adj = 0)
dev.off()

# --------------------------------
# estimate power for Healthy Homes
# --------------------------------

inputsHH <- data.frame(
  scenario = c('all_resp', 'all_cvd', 'all_mh', 'v_resp', 'v_cvd', 'v_mh'),
  outcome_rate = c(0.769, 0.455, 1.023, 3.204, 2.068, 2.842),
  outcome_disp = c(0.5, 0.5, 0.5, 5, 5, 5),
  vuln_prev = c(1, 1, 1, 0.12, 0.11, 0.18),
  study_end = 2028
)

resultsHH <- power(B = 50, inputsHH, effects = seq(0.8, 1, 0.025))

pf(inputsHH$scenario, cols = brewer.pal(6, 'Set2'), d = resultsHH)

# ------------------------------------------------
# compare to results from simple power calculation
# ------------------------------------------------

# average follow-up per household of 5.4 years (accounting for 2.18 residents per HH)

simplePower <- function (baseline = 0.769, disp = 0.5, s = 8464, effects = seq(0.7, 1, 0.001), power = 0.9) {
  p <- sapply(effects, function (x) power_NegativeBinomial(n1 = s, equal.sample = T, duration = 5.4, mu1 = baseline, mu2 = baseline * x, theta = disp)$power)
  effects[which.min(abs(p - power))]
}

inputsPASSED <- data.frame(
  s = c(8464, 8464, 8464, 1016, 931, 1523),
  baseline = c(0.769, 0.455, 1.023, 3.204, 2.068, 2.842),
  disp = c(0.5, 0.5, 0.5, 5, 5, 5) * 5.4
)

mapply(simplePower, s = inputsPASSED$s, baseline = inputsPASSED$baseline, disp = inputsPASSED$disp)
