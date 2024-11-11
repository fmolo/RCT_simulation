## -----------------------------------------------------------------------------
#| label: load_libraries
#| output: false
library(ggplot2) # for visualisations
library(scales) # for formatting axes in visualisations
library(parallel) # for parallel processing on multiple CPU cores
library(ggpubr) # for combining multiple ggplot2 graphs with shared legend
library(grid)  # for textGrob() in combined graphs


## -----------------------------------------------------------------------------
#| label: set_parameters
n <- 200 # n observations
p_max <- 1000 # p confounders
te <- 0 # treatment effect beta_t
coef_size <- 1 # beta_p (taken as constant for all confounders)
r <- c(1, 10, 10000) # number of simulations
distr <- "rnorm" # for standard-normally distributed confounders


## -----------------------------------------------------------------------------
#| label: confounder_function
# function to simulate n observations of p confounding variables of a given distribution function
generate_confounders <- function(n, p, distr = "rnorm", ...) {
    res <- sapply(1:p, FUN = function(x) {
        do.call(distr, args = list(n = n, ...))
        })
}


## -----------------------------------------------------------------------------
#| label: simulation_function
#| eval: true

simul <- function(p, r, n, te, coef_size, distr) {
    
    # initiate empty vectors to store values
    est <- numeric(r)
    ci_lower <- numeric(r)
    ci_upper <- numeric(r)
    var_between <- numeric(r)
    var_within <- numeric(r)

    for (i in 1:r) {

        # generate potential confounders
        X <- generate_confounders(n = n, p = p, distr = distr)
        
        # set effect of counfounders and random noise
        beta <- rep(coef_size, p)
        epsilon <- rnorm(n, mean = 0, sd = 1)

        # randomize treatment
        treat <- sample(rep(0:1, n / 2))

        # compute the "true" outcome y
        y <- te * treat + X %*% beta + epsilon 
      
        # fit linear regression model for treatment
        fit_lm <- lm(y ~ treat)

        # conduct analysis of variance
        fit <- anova(fit_lm)

        # record relevant values for each simulation
        est[i] <- coef(fit_lm)["treat"]
        ci_lower[i] <- confint(fit_lm)["treat", 1]
        ci_upper[i] <- confint(fit_lm)["treat", 2]
        var_between[i] <- fit["treat", "Mean Sq"]
        var_within[i] <- fit["Residuals", "Mean Sq"]
    }
    print(paste("Finished simulation for p =", p))
    print(Sys.time())
    return(list(
        est = est,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        var_between = var_between,
        var_within = var_within
    ))
}


## -----------------------------------------------------------------------------
#| label: simulate
#| eval: false

## RNGkind("L'Ecuyer-CMRG") ## needs to be set for reproducibility in parallel processing
## set.seed(0)
## 
## # using parameters set above
## # three times for three different specified numbers of simulations r
## 
## simres1 <- lapply(
##     1:p_max, simul, r = r[1], n = n, te = te,
##     coef_size = coef_size, distr = distr
## )
## 
## # parallelize the simulation across multiple CPU cores (for speed)
## # this takes long to run
## # (consider smaller p and r values to save time)
## simres2 <- mclapply(
##     1:p_max, simul, r = r[2], n = n, te = te,
##     coef_size = coef_size, distr = distr,
##     mc.cores = 6
## )
## 
## # this takes long to run
## # (consider smaller p and r values to save time)
## simres3 <- mclapply(
##     1:p_max, simul, r = r[3], n = n, te = te,
##     coef_size = coef_size, distr = distr,
##     mc.cores = 6
## )
## 
## # to reproduce the simulations without parallelization (which is probably necessary on Windows):
## # simres <- lapply(
## #     1:p_max, simul, r = r, n = n, te = te,
## #     coef_size = coef_size, distr = distr
## #     )
## 
## saveRDS(simres1, "simres_r1.rds")
## saveRDS(simres2, "simres_r2.rds")
## saveRDS(simres3, "simres_r3.rds")
## 


## -----------------------------------------------------------------------------
#| label: load_simulations
#| include: true

# alternatively we can load the results of previously run simulations from disk 
# to save time
simres1 <- readRDS("simres_r1.rds")
simres2 <- readRDS("simres_r2.rds")
simres3 <- readRDS("simres_r3.rds")


## -----------------------------------------------------------------------------
#| label: fig-figure2-panel1

plot_df1 <- data.frame(
    p = 1:p_max,
    var_between = sapply(simres1, function(x) mean(x$var_between)),
    var_within = sapply(simres1, function(x) mean(x$var_within))
    )

# Find the critical value for the ratio F for n observations and for significance level 0.05
alpha <- 0.05  # Significance level
df1 <- 1  # Between-group degrees of freedom
df2 <- n-2  # Within-group degrees of freedom
critical_f_value <- qf(1 - alpha, df1, df2)
sig_threshold <- plot_df1$var_within * critical_f_value

pointsize <- 0.8

fig2_panel1 <- ggplot(data = plot_df1, aes(x = p)) + 
        # shaded area in first layer indicates values for var_between where F > Fcrit
        geom_ribbon(aes(ymin = sig_threshold, ymax = Inf), fill = "grey80", alpha = .5) + 
        geom_point(aes(y = var_between, colour = "between groups"), size = pointsize) +
        geom_point(aes(y = var_within, colour = "within groups"), size = pointsize, alpha = 0.8) +
        #scale_y_continuous(limits = c(0, max(plot_df1$var_between))) +
        scale_colour_manual(
            name = "Variance",
            values = c(
                "between groups" = "#00608b", # blue
                "within groups" = "#ffa600" # yellow-orange
                #"between groups"= "#58508D", # purple-blue
                #"within groups" = "#FF6361" # red
                )
            ) +
        xlab("Number of potential confounders") +
        ylab("Variance") +
        theme_minimal() +
        theme(
            axis.title.x = element_text(size = 11),
            axis.text.x = element_text(size = 10, color = "black"),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.title.y = element_text(size = 10)
            ) 
        #coord_fixed(ratio = 0.025, ylim = c(0, 9000))

fig2_panel1



## -----------------------------------------------------------------------------
#| label: fig-figure2-panel2

plot_df2 <- data.frame(
    p = 1:p_max,
    var_between = sapply(simres2, function(x) mean(x$var_between)),
    var_within = sapply(simres2, function(x) mean(x$var_within))
    )

fig2_panel2 <- ggplot(data = plot_df2, aes(x = p)) + 
        geom_point(aes(y = var_between, colour = "between groups"), size = pointsize) +
        geom_point(aes(y = var_within, colour = "within groups"), size = pointsize, alpha = 0.8) +
        scale_colour_manual(
            name = "Variance",
            values = c(
                "between groups" = "#00608b", # blue
                "within groups" = "#ffa600" # yellow-orange
                #"between groups"= "#58508D", # purple-blue
                #"within groups" = "#FF6361" # red
                )
            ) +
        #scale_y_continuous(limits = c(0, max(plot_df1$var_between)/3)) +
        xlab("Number of potential confounders") +
        ylab("Variance") +
        theme_minimal() +
        theme(
            axis.title.x = element_text(size = 11),
            axis.text.x = element_text(size = 10, color = "black"),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.title.y = element_text(size = 10)
            ) 
        #coord_fixed(ratio = 0.1, ylim = c(0, 2250))

fig2_panel2


## -----------------------------------------------------------------------------
#| label: fig-figure2-panel3

plot_df3 <- data.frame(
    p = 1:p_max,
    var_between = sapply(simres3, function(x) mean(x$var_between)),
    var_within = sapply(simres3, function(x) mean(x$var_within))
    )

fig2_panel3 <- ggplot(data = plot_df3, aes(x = p)) + 
        geom_point(aes(y = var_between, colour = "between groups"), size = pointsize) +
        geom_point(aes(y = var_within, colour = "within groups"), size = pointsize/1.5) +
        scale_colour_manual(
            name = "Variance",
            values = c(
                "between groups" = "#00608b", # blue
                "within groups" = "#ffa600" # yellow-orange
                #"between groups"= "#58508D", # purple-blue
                #"within groups" = "#FF6361" # red
                )
            ) +
        #scale_y_continuous(limits = c(0, max(plot_df1$var_between)/3)) +
        xlab("Number of potential confounders") +
        ylab("Variance") +
        theme_minimal() +
        theme(
            axis.title.x = element_text(size = 11),
            axis.text.x = element_text(size = 10, color = "black"),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.title.y = element_text(size = 10)
            ) 
        #coord_fixed(ratio = 0.1, ylim = c(0, 2250))

fig2_panel3



## -----------------------------------------------------------------------------

# combine the three previously generated graphs in three panels, using ggarrange
fig2 <- ggarrange(
    fig2_panel1 + rremove("xlab") + rremove("ylab") + scale_y_continuous(
        limits = c(0, max(plot_df1$var_between)), breaks = seq(from = 0, to = 9000, by = 1000)),
    fig2_panel2 + rremove("ylab") + rremove("xlab") + scale_y_continuous(
        limits = c(0, max(plot_df1$var_between)/2.5), breaks = c(0, 1000, 2000)),  
    fig2_panel3 + rremove("xlab") + rremove("ylab") + scale_y_continuous(
        limits = c(0, max(plot_df1$var_between)/2.5), breaks = c(0, 1000, 2000)),
    ncol = 1,
    heights = c(1, 1/2.5, 1/2.5),
    common.legend = TRUE,
    legend = "top",
    labels = c("r = 1", "r = 10", "r = 10'000"),
    label.x = 0.13,
    label.y = 0.92,
    hjust = 0,
    font.label = list(family = "Helvetica", face = "bold", size = 12)
    )

# add common axis labels
fig2_annotated <- annotate_figure(
  fig2,
  left = textGrob("Variance", rot = 90, vjust = .5, hjust=0.1, gp = gpar(cex = 1, fontfamily = "Helvetica")),
  bottom = textGrob("Number of potential confounders", gp = gpar(cex= 1, fontfamily = "Helvetica"))
  )

fig2_annotated
ggsave("figures/figure2_combined.pdf", fig2_annotated, height = 8, width = 7)



## -----------------------------------------------------------------------------
expseq <- 10^seq(0, 5, by = 1)


## -----------------------------------------------------------------------------
#| eval: false
#| label: simulate_alternative

## RNGkind("L'Ecuyer-CMRG") ## needs to be set for reproducibility in parallel processing
## set.seed(1)
## 
## # this takes long to run
## simres_alt <- mclapply(
##     expseq, simul, r = r, n = n, te = te,
##     coef_size = coef_size, distr = distr,
##     mc.cores = 2
##     )
## 
## saveRDS(simres_alt, "simres_alt.rds")


## -----------------------------------------------------------------------------
#| eval: true
#| label: load_simulations_alternative

# alternatively we can load the results of previously run simulations from disk
simres_alt <- readRDS("simres_alt.rds")


## -----------------------------------------------------------------------------
#| label: extract_metrics

# extract empirical bias and variance of the estimator
bias <- sapply(simres_alt, function(x) mean(x$est) - te)
variance <- sapply(simres_alt, function(x) var(x$est))

# bias and variance adjusted for Monte Carlo error
bias_mc_error <- sapply(simres_alt, function(x) sd(x$est) / sqrt(r))
variance_mc_error <- sapply(simres_alt, function(x) sd(x$est)^2 * sqrt(2 / (r - 1)))

# calculate coverage of the confidence interval
get_coverage <- function(ci_lower, ci_upper, te) {
    sum(ci_lower < te & ci_upper > te) / length(ci_lower)
}

coverage <- sapply(
    simres_alt, function(x) get_coverage(x$ci_lower, x$ci_upper, te)
    )
coverage_mc_error <- sqrt(coverage * (1 - coverage) / r)

plot_df_alt <- data.frame(
    p = expseq,
    bias = bias,
    bias_mc_error = bias_mc_error,
    variance = variance,
    variance_mc_error = variance_mc_error,
    coverage = coverage,
    coverage_mc_error = coverage_mc_error
)


## -----------------------------------------------------------------------------
#| eval: true
#| label: fig-bias
#| fig.cap: Empirical bias of the estimator with growing number of confounders. The x-axis is on a logarithmic scale. The dot represents the average bias across simulations, the bar spans +/-1 Monte Carlo standard error (i.e. the standard deviation across simulations, enlarged to include the imprecision due to the finite number of simulations.)

options(scipen = 999) # force ggplot to not use scientific notation for the x-axis

ggplot(plot_df_alt, aes(x = p, y = bias)) +
    geom_pointrange(aes(
        ymin = bias - bias_mc_error,
        ymax = bias + bias_mc_error
        ), fatten = 1) +    
    scale_x_log10(breaks = log_breaks(n = 6), ) +
    theme_minimal()



## -----------------------------------------------------------------------------
#| eval: true
#| label: fig-variance
#| fig.cap: Variance of the estimator with growing number of confounders. The x-axis is on a logarithmic scale. The dot represents the average variance across simulations, the bar spans +/-1 Monte Carlo standard error.
ggplot(plot_df_alt, aes(x = p, y = variance)) +
    geom_pointrange(aes(
        ymin = variance - variance_mc_error,
        ymax = variance + variance_mc_error
        ), fatten = 1) +
    scale_x_log10(breaks = log_breaks(n = 6)) +
    theme_minimal()


## -----------------------------------------------------------------------------
#| eval: true
#| label: fig-coverage
#| fig.cap: Coverage of the 95% confidence interval with growing number of confounders. The x-axis is on a logarithmic scale. The dot represents the average coverage across simulations. The bar spans +/-1 Monte Carlo standard error, but this is so small it is barely visible.

ggplot(plot_df_alt, aes(x = p, y = coverage)) +
    geom_pointrange(aes(
        ymin = coverage - coverage_mc_error,
        ymax = coverage + coverage_mc_error
        ), fatten = 1) +
    scale_x_log10(breaks = log_breaks(n = 6)) +
    scale_y_continuous(
        limits = c(0, 1),
        labels = scales::percent,
        breaks = c(0, 0.25, 0.5, 0.75, 0.95, 1)
        ) +
    theme_minimal() +
    theme(panel.grid.minor.y = element_blank())



## -----------------------------------------------------------------------------
#| label: session_info
sessionInfo()

