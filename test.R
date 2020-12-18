N = 1000
Z = rnorm(N, 0, 1)
Y = rnorm(N, 0, 1) / 2 + Z / 2
X = rnorm(N, 0, 1) / 2 + Z / 2

mod = lm(Y ~ X + Z)

r = resid(mod)

cor(r, Z)

# simplest model: x correlated with y and also with z

# the "network" here is a process for pairing up predictions / observations


# Y is fixed; noisy measurement X; noise correlated with network
# intuition: if you learn a model Y ~ X, and the residuals are network
# correlated, then you will get a biased homophily measurement
# EVEN THOUGH you get a unbiased individual level prediction
# 
# The way in which our signal is noisy is network-correlated
# e.g. the name example
# 
# Map: draw Y. Draw X = f(Y, X_tilde) where X_tilde is noise.
# Generate the network according to both Y and X_tilde
# Then, even though the predictions are good among *nodes* they
# will be poor among *edges* for the purpose of quantification

N = 1000
Y = rnorm(N, 0, 1)
X = Y + rnorm(N, 0, 1)

mod = lm(Y ~ X)