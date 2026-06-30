# ==============================================================================
# Exploratory Script: 5-Point SIG Optimal Design for Norberg Kinetics (Temperature)
# ==============================================================================
# This script is designed for interactive exploration. It calculates a 5-point 
# optimal experimental design using the acebayes package, plots the correlated 
# parameter priors, and visualizes the resulting design. We strongly recommend 
# reading the overview of ACE in the Supplementary Information of the paper.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Setup and Package Loading
# ------------------------------------------------------------------------------
# We load acebayes for the design calculation, MASS/Matrix for correlated priors,
# and GGally for the matrix pair plots.

# Check for required packages and install if missing
local({
  pkgs <- c(
    "acebayes",
    "dplyr",
    "ggplot2",
    "here",
    "parallel",
    "tidyr",
    "GGally",
    "MASS",
    "Matrix"
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
library(GGally)
library(MASS)
library(Matrix)

# This is a safety check to make sure that this script is 
# in the associated project folder.
here::i_am("temperature_SIG_design_exploration.R")

# Source the custom C++ functions for stability
source(here("00_support_file_stable_SIG_CPP_function.R"))
source(here("00_support_file_stable_utility_function.R"))
environment(utilitynlmTemp) <- asNamespace("acebayes")
assignInNamespace("utilitynlm", utilitynlmTemp, ns = "acebayes")

# Source the temperature prior file
# We are using a custom function for temperature priors here
# The function enforces correlations between some parameters
# And gives us a joint probability distribution of the parameters
# If the parameters are independent, they often result in unrealistic TPCs
# The script with the temperature prior function is in the Github repo.
# Feel free to customize it and calculate new designs.
source(here("00_support_file_temperature_prior.R"))

# The function inside the sourced file is named 'prior', we rename it 
# to 'norberg_prior'
norberg_prior <- prior 

# Set seed for reproducibility of the starting designs during exploration
set.seed(42)

# ------------------------------------------------------------------------------
# 2. Visualize Correlated Prior Distribution 
# ------------------------------------------------------------------------------
# The Norberg model parameters (a, b, tmax, tmin) are correlated. We use 
# GGally::ggpairs to plot the distributions and their correlations.

# Draw 2,000 samples (reduced from 10k so the scatterplot matrix renders quickly)
prior_samples <- as.data.frame(norberg_prior(2000))

# --- Plot the Priors Matrix ---
cat("\nRendering correlation matrix plot for priors (this may take a few seconds)...\n")

p_matrix <- ggpairs(
  prior_samples,
  lower = list(continuous = wrap("points", alpha = 0.15, size = 0.5, color = "#0072B2")),
  diag = list(continuous = wrap("densityDiag", fill = "#56B4E9", alpha = 0.6)),
  upper = list(continuous = wrap("cor", size = 4, color = "black")),
  title = "Prior distributions and correlations for Norberg parameters"
) + 
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 10))

print(p_matrix)

# ------------------------------------------------------------------------------
# 3. Setup Design Space and ACE Algorithm parameters
# ------------------------------------------------------------------------------
# Range of temperatures at which an experimental unit can be placed
low <- 5
upp <- 35

# Number of experimental units in the design
n   <- 5   # We are strictly looking for a 5-point design

# Grid of possible temperatures: 5 to 35 with length.out = 61 (step size of 0.5 deg C)
grid_pts <- seq(from = low, to = upp, length.out = 61)

limits <- function(d, i, j){
  grid_pts
}

# ------------------------------------------------------------------------------
# 4. Calculate Optimal Design
# ------------------------------------------------------------------------------
# First we try a single starting design and look at the final design we get
start.d <- matrix(sample(grid_pts, n),
                  nrow = n, ncol = 1,
                  dimnames = list(as.character(1:n), c("temp")))

# acenlm finds the SIG-optimal design (criterion = "SIG")
# (N1, N2, B) are reduced for rapid exploration as compared to the manuscript.
temp_5_single_start <- acenlm(
  formula = ~ (a*exp(b*temp))*(tmax - temp)*(temp - tmin),
  start.d = start.d, 
  prior   = norberg_prior,
  N1 = 15, N2 = 50, B = c(5000, 1000), # Reduced from 40, 200, c(40000, 2000)
  lower = low, upper = upp, limits = limits,
  method = "MC", criterion = "SIG", progress = TRUE
)

# Look at the summary
temp_5_single_start

# Look at whether the algorithm converged in Phase 1
plot(temp_5_single_start$phase1.trace, 
     type = "b", pch = 19, col = "blue",         
     xlab = "Iteration (Phase I)", ylab = "Approximate Expected SIG",
     main = "ACE Algorithm Convergence Trace (Phase I)")

####### Mac/Linux version: Designs from several starting designs ########
# If you are on a Windows machine skip ahead to the Windows block

# For exploration, we start with only a few starting designs (C) 
C <- 4 
start.d <- list()
for(i in 1:C){
  start.d[[i]] <- matrix(sample(grid_pts, n),
                         nrow = n, ncol = 1,
                         dimnames = list(as.character(1:n), c("temp")))
}

cores <- max(1, detectCores() - 1)

temp_5_multiple_starts <- pacenlm(
  formula = ~ (a*exp(b*temp))*(tmax - temp)*(temp - tmin),
  start.d = start.d, 
  prior   = norberg_prior,
  N1 = 15, N2 = 50, B = c(5000, 1000),
  lower = low, upper = upp, limits = limits,
  method = "MC", criterion = "SIG", mc.cores = cores
)

# This is the best final design.
temp_5_multiple_starts$d

# These are ALL the final designs.
temp_5_multiple_starts$final.d

# For the best final design, look at whether Phase I converged
plot(temp_5_multiple_starts$phase1.trace, 
     type = "b",           
     pch = 19,             
     col = "blue",         
     xlab = "Iteration (Phase I Stage)", 
     ylab = "Approximate Expected SIG",
     main = "ACE Algorithm Convergence Trace")

# For the best final design, look at whether Phase II converged
plot(toxin_5_multiple_starts$phase2.trace, 
     type = "b",           
     pch = 19,             
     col = "blue",         
     xlab = "Iteration (Phase I Stage)", 
     ylab = "Approximate Expected SIG",
     main = "ACE Algorithm Convergence Trace")


####### Windows version: Designs from several starting designs ########
# Mac/Linux users: use the pacenlm block above instead.
# Windows users: comment out the pacenlm block and use this.
# ==============================================================================

# cores <- max(1, detectCores() - 1)
# cl <- makeCluster(cores)
# 
# # Load acebayes and MASS/Matrix on all workers
# clusterEvalQ(cl, {
#   library(acebayes)
#   library(MASS)
#   library(Matrix)
# })
# 
# clusterExport(cl,
#               varlist = c("norberg_prior", "low", "upp", "n", "limits", "grid_pts"),
#               envir   = environment())
# 
# # Source the stable versions of utility calculations on workers.
# clusterEvalQ(cl, {
#   source(here::here("support_file_stable_SIG_CPP_function.R"))
#   source(here::here("support_file_stable_utility_function.R"))
#   environment(utilitynlmTemp) <- asNamespace("acebayes")
#   assignInNamespace("utilitynlm", utilitynlmTemp, ns = "acebayes")
# })
# 
# run_one_start <- function(sd) {
#   acenlm(
#     formula = ~ (a*exp(b*temp))*(tmax - temp)*(temp - tmin),
#     start.d = sd,
#     prior   = norberg_prior,
#     N1 = 15, N2 = 50, B = c(5000, 1000),
#     lower = low, upper = upp, limits = limits,
#     method = "MC", criterion = "SIG", progress = FALSE
#   )
# }
# 
# ace_runs <- parLapply(cl, start.d, run_one_start)
# 
# n.assess <- 20
# B_assess <- 5000   
# 
# assess_one <- function(run, n.assess, B_assess) {
#   u <- run$utility                 
#   d <- run$phase2.d                
#   ev <- numeric(n.assess)
#   for (k in seq_len(n.assess)) ev[k] <- mean(u(d = d, B = B_assess))
#   ev
# }
# 
# clusterExport(cl, varlist = c("assess_one", "n.assess", "B_assess"),
#               envir = environment())
# eval_list <- parLapply(cl, ace_runs, assess_one, n.assess, B_assess)
# 
# stopCluster(cl)   
# 
# besti <- which.max(vapply(eval_list, mean, numeric(1)))
# 
# temp_5_multiple_starts <- list(
#   d            = ace_runs[[besti]]$phase2.d,                 
#   final.d      = lapply(ace_runs, function(r) r$phase2.d),   
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
SIG_opt_design <- temp_5_multiple_starts$d

design_df <- data.frame(S = as.numeric(SIG_opt_design), Design = "Optimal (n=5)")

designs_agg <- design_df %>% 
  group_by(S, Design) %>% 
  summarise(count = n(), .groups = 'drop')

stack_step <- 0.025
y_base <- 0

dot_df <- designs_agg %>%
  uncount(count, .id = "stack_id") %>%
  mutate(y = y_base + (stack_id - 1) * stack_step)

# --- Plot the Design ---
p_design <- ggplot() +
  geom_hline(yintercept = -0.05, color = "black", linewidth = 0.5) +
  geom_point(data = dot_df, aes(x = S, y = y),
             fill = "#0072B2", shape = 21, alpha = 0.85, 
             size = 6, colour = "white", stroke = 0.5) +
  scale_x_continuous(limits = c(0, 40), breaks = seq(0, 40, 5)) +
  coord_cartesian(ylim = c(-0.1, max(dot_df$y) + 0.2)) +
  labs(title = "5-point SIG-optimal design for Temperature (Norberg)",
       subtitle = "Each point represents one experimental unit",
       x = expression(Temperature~(degree*C)), 
       y = NULL) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(size = 16),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

print(p_design)
