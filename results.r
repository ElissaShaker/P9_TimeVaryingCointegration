library(ggplot2)
library(data.table)

# dim of chebyshev polynomial
cat("Chosen m (by Hannan-Quinn):", m, "\n")
cat("k (number of series):", k, "   p (VAR lag):", p, "   n-p-1:", n-p-1, "\n\n")

# Show lrt stats and p-values for all r (1..k) at the chosen m
for(r in 1:k){
  stat <- lrtvc[m, r]
  pval <- lrtvcpv[m, r]
  df <- m * r * k
  decision <- ifelse(is.na(pval), "NA", ifelse(pval < 0.05, "Reject H0 (time-varying)", "Fail to reject H0 (constant)"))
  cat(sprintf("r = %d : LR stat = %.4f  df = %d  p-value = %.4g   => %s\n", r, stat, df, pval, decision))
}

# Optional: print full matrices for inspection
print("statistic (first few rows):")
print(lrtvc[1:min(10,nrow(lrtvc)), , drop = FALSE])
print("p-value (first few rows):")
print(lrtvcpv[1:min(10,nrow(lrtvcpv)), , drop = FALSE])


################################
# Assume betat is your n_obs x (k*m) matrix from tvcoint
# k = 2 (land, sea), m = chosen Chebyshev dimension (here 3)
n_obs <- nrow(betat)
time_index <- seq(1, n_obs)  # or use actual dates if available

# For r = 1, extract the first k columns of betat for the long-run relationship
beta_r1 <- betat[, 1:k]

# Convert to long format for ggplot
beta_dt <- data.table(
  time = rep(time_index, each = k),
  series = rep(c("Land", "Sea"), times = n_obs),
  beta = as.vector(beta_r1)
)

# Plot
ggplot(beta_dt, aes(x = time, y = beta, color = series)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Time-Varying Cointegrating Vector (r=1)",
    x = "Time Index",
    y = expression(beta[t]),
    color = "Series"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  )

