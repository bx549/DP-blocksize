## implementation of value iteration for example 1.3.1 in
## bertsekas, dynamic programming and optimal control volume 2.
## this code is for illustrative purposes. it can serve as a
## starting point to solve other DP problems.
library(tidyverse)

S <- 1:2  # state space
C <- 1:2  # controls

## transition probabilities
P <- array(data = c(3/4, 3/4, 1/4, 1/4,
                  1/4, 1/4, 3/4, 3/4),
           dim = c(2, 2, 2),
           dimnames = list(statei=as.character(S),
                            statej=as.character(S),
                            control=as.character(C)))

## cost to apply control u in state i
g <- array(data = c(2, 1, 0.5, 3),
           dim = c(2, 2),
           dimnames = list(state=as.character(S),
                           control=as.character(C)))

alpha <- 0.9   # discount factor

## applies the DP iteration for state i
FJ <- function(i) {
    cost <- rep(Inf, length(C))   # holds the cost for each control
    for (u in C) {
        cost.to.go <- 0   # holds expected cost-to-go for control u
        for (j in S) {
            cost.to.go <- cost.to.go + P[i,j,u]*J[j]
        }
        cost[u] <- g[i,u] + alpha*cost.to.go
    }
    cost.star <- min(cost)
    ctrl <- which(near(cost.star, cost))
    ## return the min cost and the control that acheived
    ## the min cost
    list(cost.star, ctrl)
}

## apply the value iteration algorithm.
## note that this is the gauss-seidel version of value iteration
## because we iterate one state at a time and use the interim results,
## hence the notation FJ(i) as in bertsekas.
eps <- .0001           # tolerance for convergence
k <- 0                 # iteration number
J <- rep(0, length(S)) # initially, start the value iteration with zero cost
J.prev <- rep(eps + 1, length(S)) # to check convergence
policy <- rep(NA, length(S))

while (!all(near(J - J.prev, 0, tol=eps))) {
    k <- k + 1
    J.prev <- J
    for (i in 1:length(S)) {
        lst <- FJ(i)
        J[i] <- lst[[1]]
        policy[i] <- lst[[2]]
    }
}
