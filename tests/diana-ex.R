library(cluster)
options(digits = 6)
data(votes.repub)
summary(diana(votes.repub, metric = "manhattan", stand = TRUE))
summary(diana(daisy(votes.repub), diss = TRUE))

data(agriculture)
summary(diana(agriculture))

data(ruspini)
summary(diana(ruspini))
summary(diana(ruspini, metric = "manhattan"))
