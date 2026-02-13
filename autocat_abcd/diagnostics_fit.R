# 1. Load the data
df <- read.csv("protocell_sim_prueba.csv")

n_total <- nrow(df)
sixtysix <- floor(n_total / 2)
df_short <- df[sixtysix:n_total, ]

# 2. Fit the linear model (Growth Rate alpha)
# We model log_volume as a function of time
myAdditive <- lm(log_volume ~ time, data = df_short)

# 3. View the summary (Alpha, R^2, and p-values)
summary(myAdditive)

# 4. The "Magic" Diagnostics
# This creates the 2x2 grid you remember
oldpar <- par(oma=c(0,0,3,0), mfrow=c(2,2))

plot(myAdditive)

# Add a main title to the outer margin

par(oldpar)

library(moments)
kurtosis(residuals(myAdditive))


plot(volume ~ time, data= df)
