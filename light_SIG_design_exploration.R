# ==============================================================================
# Exploratory Script: 5-Point SIG Optimal Design for Eilers-Peeters Kinetics (Light)
# ==============================================================================
# This script is designed for interactive exploration. It calculates a 5-point 
# optimal experimental design using the acebayes package, plots the parameter 
# priors and visualizes the resulting design. We strongly recommend reading the
# overview of ACE in the Supplementary Information of the paper before
# exploring designs with this script. 
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Setup and Package Loading
# ------------------------------------------------------------------------------
# We use the R package acebayes to do the design calculation
# We load acebayes for the design calculation, and the tidyverse/patchwork 
# suite for design plotting.

# Check for required packages and install if missing
local({
  pkgs <- c(
    "acebayes",
    "dplyr",
    "ggplot2",
    "here",
    "parallel",
    "patchwork",
    "tidyr"
  )
  
  missing <- pkgs[
    !vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
  ]
  
  if (length(missing)) {
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
})

# Load libraries
library(here)
library(acebayes)
library(parallel)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# This is a safety check to make sure that this script is 
# in the associated project folder. 
here::i_am("light_SIG_design_exploration.R")

# The original package's code occasionally results in Inf or -Inf SIG values 
# This happens due to numerical issues: underflow in case of 
# large negative log-likelihood values.
# We have fixed this issue and are sourcing the modified code below
# The files are in the Github repo. Save them into the same folder that you save this file.
source(here("00_support_file_stable_SIG_CPP_function.R"))
source(here("00_support_file_stable_utility_function.R"))
environment(utilitynlmTemp) <- asNamespace("acebayes")
assignInNamespace("utilitynlm", utilitynlmTemp, ns = "acebayes")

# Set seed for reproducibility of the starting designs during exploration
# Note that acebayes does not allow for reproduction of the final designs
set.seed(42)

# ------------------------------------------------------------------------------
# 2. Define Prior Distribution & Visualize
# ------------------------------------------------------------------------------
# The Eilers-Peeters model parameters (mu_max, alpha, i_opt, and measurement error 
# variance sig2) are drawn from lognormal distributions.

eilers_peeters_prior <- function(B){
  mu_max <- rlnorm(n = B, meanlog = 0.01, sdlog = 0.4)
  alpha  <- rlnorm(n = B, meanlog = -3, sdlog = 0.8)
  i_opt  <- rlnorm(n = B, meanlog = 5.5, sdlog = 0.3)
  sig2   <- rlnorm(n = B, meanlog = -2.3*2, sdlog = 0.1*2)
  
  out <- cbind(mu_max, alpha, i_opt, sig2)
  colnames(out) <- c("mu_max", "alpha", "i_opt", "sig2")
  return(out)
}

# --- Plot the Priors ---
# Draw 10,000 samples to visualize what the prior space looks like
prior_samples <- as.data.frame(eilers_peeters_prior(10000))

p_mu <- ggplot(prior_samples, aes(x = mu_max)) +
  geom_density(fill = "#56B4E9", alpha = 0.6, color = NA) +
  coord_cartesian(xlim = c(0, 5)) +
  labs(title = expression("Prior for" ~ mu[max]), x = expression(mu[max]), y = "Density") +
  theme_minimal(base_size = 14)

p_alpha <- ggplot(prior_samples, aes(x = alpha)) +
  geom_density(fill = "#009E73", alpha = 0.6, color = NA) +
  coord_cartesian(xlim = c(0, 0.3)) +
  labs(title = expression("Prior for" ~ alpha), x = expression(alpha), y = "Density") +
  theme_minimal(base_size = 14)

p_iopt <- ggplot(prior_samples, aes(x = i_opt)) +
  geom_density(fill = "#D55E00", alpha = 0.6, color = NA) +
  coord_cartesian(xlim = c(0, 800)) +
  labs(title = expression("Prior for" ~ I[opt]), x = expression(I[opt]), y = "Density") +
  theme_minimal(base_size = 14)

p_sig2 <- ggplot(prior_samples, aes(x = sig2)) +
  geom_density(fill = "#CC79A7", alpha = 0.6, color = NA) +
  coord_cartesian(xlim = c(0, 0.03)) +
  labs(title = expression("Prior for" ~ sigma^2), x = expression(sigma^2), y = "Density") +
  theme_minimal(base_size = 14)

# Show prior plots in a 2x2 grid using patchwork
print((p_mu | p_alpha) / (p_iopt | p_sig2))

# ------------------------------------------------------------------------------
# 3. Setup Design Space and ACE Algorithm parameters
# ------------------------------------------------------------------------------
# Range of light intensities at which an experimental unit can be placed
low <- 10
upp <- 1000

# Number of experimental units in the design
n   <- 5   # We are strictly looking for a 5-point design

# Experimental units can't typically be placed at any arbitrary value
# Limits function to constrain the acebayes search space to a grid.
# This grid has 100 evenly spaced points from 10 to 1000.
limits <- function(d, i, j){
  seq(from = low, to = upp, length.out = 100)
}

# ------------------------------------------------------------------------------
# 4. Calculate Optimal Design
# ------------------------------------------------------------------------------
# First we try a single starting design and look at the final design we get
# Define a starting design by randomly picking 5 points from the grid
start.d <- matrix(sample(seq(from = low, to = upp, length.out = 100), n),
                  nrow = n, ncol = 1,
                  dimnames = list(as.character(1:n), c("light")))

# acenlm finds the SIG-optimal design (criterion = "SIG"),
# averaging over priors using Monte Carlo (method = "MC").
# (N1, N2, B) are reduced for rapid exploration as compared to the manuscript.
light_5_single_start <- acenlm(
  formula = ~ ((mu_max * light) / (((mu_max / (alpha * i_opt^2)) * (light^2)) + 
                                     ((1 - ((2 * mu_max) / (alpha * i_opt))) * light) + 
                                     (mu_max / alpha))),
  start.d = start.d, 
  prior   = eilers_peeters_prior,
  N1 = 15, N2 = 50, B = c(5000, 1000), # Reduced from 40, 200, c(40000, 2000)
  lower = low, upper = upp, limits = limits,
  method = "MC", criterion = "SIG", progress = TRUE
)


# Look at the summary
light_5_single_start

# Look at the design at the end of Phase 1
light_5_single_start$phase1.d

# Look at whether the algorithm converged
plot(light_5_single_start$phase1.trace, 
     type = "b",           # "b" means plot BOTH points and lines
     pch = 19,             # pch = 19 gives solid circle points
     col = "blue",         # Line and point color
     xlab = "Iteration (Phase I)", 
     ylab = "Approximate Expected SIG",
     main = "ACE Algorithm Convergence Trace")

# Look at the design at the end of Phase 2 (final design)
light_5_single_start$phase2.d

# Look at whether Phase II converged
plot(light_5_single_start$phase2.trace, 
     type = "b",           
     pch = 19,             
     col = "blue",         
     xlab = "Iteration (Phase II)", 
     ylab = "Approximate Expected SIG",
     main = "ACE Algorithm Convergence Trace")

####### Mac/Linux version: Designs from several starting designs ########
# If you are on a Windows machine skip ahead to the Windows block

# For exploration, we start with only a few starting designs (C) 
C <- 4 
start.d <- list()
for(i in 1:C){
  start.d[[i]] <- matrix(sample(seq(from = low, to = upp, length.out = 100), n),
                         nrow = n, ncol = 1,
                         dimnames = list(as.character(1:n), c("light")))
}

# Number of cores to use (leaves 1 free to keep the computer usable)
cores <- max(1, detectCores() - 1)

# As with acenlm above, (N1, N2, B) are reduced for rapid exploration.
light_5_multiple_starts <- pacenlm(
  formula = ~ ((mu_max * light) / (((mu_max / (alpha * i_opt^2)) * (light^2)) + 
                                     ((1 - ((2 * mu_max) / (alpha * i_opt))) * light) + 
                                     (mu_max / alpha))),
  start.d = start.d, 
  prior   = eilers_peeters_prior,
  N1 = 15, N2 = 50, B = c(5000, 1000),
  lower = low, upper = upp, limits = limits,
  method = "MC", criterion = "SIG", mc.cores = cores
)

# Look at the summary, including how long it took.
light_5_multiple_starts

# This is the best final design.
light_5_multiple_starts$d

# These are ALL the final designs.
light_5_multiple_starts$final.d


# For the best final design, look at whether Phase I converged
plot(light_5_multiple_starts$phase1.trace, 
     type = "b",           
     pch = 19,             
     col = "blue",         
     xlab = "Iteration (Phase I Stage)", 
     ylab = "Approximate Expected SIG",
     main = "ACE Algorithm Convergence Trace")

# For the best final design, look at whether Phase II converged
plot(light_5_multiple_starts$phase2.trace, 
     type = "b",           
     pch = 19,             
     col = "blue",         
     xlab = "Iteration (Phase I Stage)", 
     ylab = "Approximate Expected SIG",
     main = "ACE Algorithm Convergence Trace")

####### Windows version: Designs from several starting designs ########
# pacenlm uses mclapply (fork-based), which does not work on Windows.
# This block reproduces pacenlm by running acenlm on each starting design
# across a PSOCK cluster (parLapply), then selecting the best final design
# the same way pacenlm does: re-evaluate each design n.assess times and
# compare the MEAN estimated SIG (each SIG value is a noisy MC estimate).
#
# Mac/Linux users: use the pacenlm block above instead.
# Windows users: comment out the pacenlm block and use this.
# ==============================================================================

# cores <- max(1, detectCores() - 1)
# 
# # 1. Create a PSOCK cluster (works on Windows, Mac, and Linux)
# cl <- makeCluster(cores)
# 
# # 2. PSOCK workers start empty: load packages and ship needed objects to each.
# #    (Forking would inherit these automatically; PSOCK does not.)
# clusterEvalQ(cl, library(acebayes))
# clusterExport(cl,
#               varlist = c("eilers_peeters_prior", "low", "upp", "n", "limits"),
#               envir   = environment())
# 
# # 3. Source the stable versions of utility calculations.
# clusterEvalQ(cl, {
#   source(here::here("1_scripts", "support_file_stable_SIG_CPP_function.R"))
#   source(here::here("1_scripts", "support_file_stable_utility_function.R"))
#   environment(utilitynlmTemp) <- asNamespace("acebayes")
#   assignInNamespace("utilitynlm", utilitynlmTemp, ns = "acebayes")
# })
# 
# 
# # 4. Run acenlm once per starting design, in parallel across the cluster.
# run_one_start <- function(sd) {
#   acenlm(
#     formula = ~ ((mu_max * light) / (((mu_max / (alpha * i_opt^2)) * (light^2)) + 
#                                        ((1 - ((2 * mu_max) / (alpha * i_opt))) * light) + 
#                                        (mu_max / alpha))),
#     start.d = sd,
#     prior   = eilers_peeters_prior,
#     N1 = 15, N2 = 50, B = c(5000, 1000),
#     lower = low, upper = upp, limits = limits,
#     method = "MC", criterion = "SIG", progress = FALSE
#   )
# }
# 
# ace_runs <- parLapply(cl, start.d, run_one_start)
# 
# # 5. Select the best final design the way pacenlm does.
# n.assess <- 20
# B_assess <- 5000   # = B[1] used above
# 
# assess_one <- function(run, n.assess, B_assess) {
#   u <- run$utility                 # the design-scale utility acenlm returns
#   d <- run$phase2.d                # the final (Phase II) design for this start
#   ev <- numeric(n.assess)
#   for (k in seq_len(n.assess)) ev[k] <- mean(u(d = d, B = B_assess))
#   ev
# }
# 
# # Ship the assessment helper + args, then evaluate in parallel.
# clusterExport(cl, varlist = c("assess_one", "n.assess", "B_assess"),
#               envir = environment())
# eval_list <- parLapply(cl, ace_runs, assess_one, n.assess, B_assess)
# 
# stopCluster(cl)   # always shut the cluster down when done
# 
# # 7. Best start = highest mean estimated SIG, and rebuild a pacenlm-like object
# besti <- which.max(vapply(eval_list, mean, numeric(1)))
# 
# #This returns an equivalent object back
# light_5_multiple_starts <- list(
#   d            = ace_runs[[besti]]$phase2.d,                 # best final design
#   final.d      = lapply(ace_runs, function(r) r$phase2.d),   # all final designs
#   phase1.d     = ace_runs[[besti]]$phase1.d,
#   phase2.d     = ace_runs[[besti]]$phase2.d,
#   phase1.trace = ace_runs[[besti]]$phase1.trace,
#   phase2.trace = ace_runs[[besti]]$phase2.trace,
#   eval         = eval_list[[besti]],
#   besti        = besti
# )
# 
# ------------------------------------------------------------------------------
# 5. Extract and Visualize the Resulting Design
# ------------------------------------------------------------------------------
# Extract the best design points from the acebayes object
SIG_opt_design <- light_5_multiple_starts$d

# Prepare design to be plotted
design_df <- data.frame(S = as.numeric(SIG_opt_design), Design = "Optimal (n=5)")

# Aggregate counts for each unique level
designs_agg <- design_df %>% 
  group_by(S, Design) %>% 
  summarise(count = n(), .groups = 'drop')

# Uncount and map to a Y-axis stack position
# Adjust the vertical distance between two points stacked on top of each other
stack_step <- 0.025
y_base <- 0

dot_df <- designs_agg %>%
  uncount(count, .id = "stack_id") %>%
  mutate(y = y_base + (stack_id - 1) * stack_step)

# --- Plot the Design ---
p_design <- ggplot() +
  # Draw a baseline to anchor the points
  geom_hline(yintercept = -0.05, color = "black", linewidth = 0.5) +
  geom_point(data = dot_df, aes(x = S, y = y),
             fill = "#0072B2", shape = 21, alpha = 0.85, 
             size = 6, colour = "white", stroke = 0.5) +
  scale_x_continuous(limits = c(0, 1050), breaks = seq(0, 1000, 200)) +
  coord_cartesian(ylim = c(-0.1, max(dot_df$y) + 0.2)) +
  labs(title = "5-point SIG-optimal design for light (Eilers-Peeters)",
       subtitle = "Each point represents one experimental unit",
       x = expression(Light~intensity~(mu*mol~m^{-2}~s^{-1})), 
       y = NULL) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(size = 16),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

# Display the final design plot
print(p_design)
