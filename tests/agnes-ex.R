library(cluster)
options(digits = 6)
data(votes.repub)
agn1 <- agnes(votes.repub, metric = "manhattan", stand = TRUE)
summary(agn1)
agn2 <- agnes(daisy(votes.repub), diss = TRUE, method = "complete")
summary(agn2)

data(agriculture)
summary(agnes(agriculture))

data(ruspini)
summary(agnes(ruspini))
summary(agnes(ruspini, metric = "manhattan"))

