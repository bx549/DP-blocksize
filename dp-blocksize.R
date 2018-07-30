## value iteration for dynamic block size
library(tidyverse)

## state space is the number of txs in the mempool
lb <- 1000
ub <- 10000
S <- seq(lb, ub, by=1)
## control is the block size in MB
## 1 MB gets 2000 txs
## 2 MB gets 4000 txs
C <- 1:2

lambda <- 2500  # arrival rate of new txs (per 10 minutes)
c <- .10        # user delay cost in $ per minute

plot(S, dpois(S, 2500))

f <- function(x, u, w) {   # system evolution
    stopifnot(u %in% 1:2)
    x - u*2000 + w
}

## transition probabilities
P <- array(0,
           dim = c(length(S), length(S), length(C)),
           dimnames = list(statei=as.character(S),
                           statej=as.character(S),
                           control=as.character(C)))
for (i in 1:length(S)) {
#    for (j in 1:length(S)) {
        for (u in 1:length(C)) {
            ## compute transition prob. in order to go from i to j
            ## when u is applied, w would need to be j-i+u*2000
            w <- pmax(S - S[i] + C[u]*2000, 0)
            P[i,,u] <- dpois(w, lambda)
        }
#    }
}

## cost to apply control u in state i
g <- array(data = 0,
           dim = c(length(S), length(C)),
           dimnames = list(state=as.character(S),
                           control=as.character(C)))

for (i in 1:length(S)) {
    for (u in 1:length(C)) {
        g[i,u] <- max(S[i] - C[u]*2000*c*10, 0)
    }
}

## applies the DP iteration for state i. note that the argument i is
## is the index into the state space S.
## note: the error bounds are not yet implemented
FJ <- function(i) {
    cost <- rep(Inf, length(C))   # holds the cost for each control
    for (u in 1:length(C)) {
        cost.to.go <- sum( P[i,,u]*J )
        cost[u] <- g[i,u] + alpha*cost.to.go
    }
    cost.star <- min(cost)
    ctrl <- which(near(cost.star, cost)) # store the ctrl that achieves the min cost
    ## note that there could be a tie and length(ctrl) > 1, hence just choose 1st entry
    list(cost.star, ctrl[1])
}

## apply the value iteration algorithm.
eps <- .0001                      # tolerance for convergence
k <- 0                            # iteration number
J <- rep(100, length(S))          # initial starting values
J.prev <- rep(eps + 1, length(S)) # to check convergence
policy <- rep(NA, length(S))

while (!all(near(J - J.prev, 0, tol=eps))) {
    k <- k + 1
    J.prev <- J
    for (j in 1:length(S)) {
        lst <- FJ(j)    # pass the index
        J[j] <- lst[[1]]
        policy[j] <- lst[[2]]
        if (near(k %% 100, 0)) { # print a message to check progress
            message("k = ", k, " J(", j, ") = ", J[j])
        }
    }
}
