## implementation of value iteration for example 1.3.1 in
## bertsekas, dynamic programming and optimal control volume 2.
library(tidyverse)

S <- 1:2
C <- 1:2

P <- array(data = c(3/4, 3/4, 1/4, 1/4,
                  1/4, 1/4, 3/4, 3/4),
           dim = c(2, 2, 2),
           dimnames = list(statei=as.character(S),
                            statej=as.character(S),
                            control=as.character(C)))

g <- array(data = c(2, 1, 0.5, 3),
           dim = c(2, 2),
           dimnames = list(state=as.character(S),
                           control=as.character(C)))

alpha <- 0.9

## applies the DP iteration for state i
## note: the error bounds are not yet implemented
FJ <- function(i) {
    cost <- rep(Inf, length(C))   # holds the cost for each control
    for (u in C) {
        cost.to.go <- 0   # holds expected cost-to-go for control u
        for (j in S) {
            cost.to.go <- cost.to.go + P[i,j,u]*J[j]
        }
        cost[u] <- g[i,u] + alpha*cost.to.go
    }
    min(cost)
}


## apply the value iteration algorithm.
## note that this is the gauss-seidel version of value iteration
## because we iterate one state at a time and use the interim results,
## hence the notation FJ(i).
## if we implment the error bounds then we would need to
## use the ordinary version of value iteration where we iterate
## for all states simultaneously, that is TJ(i) if we follow
## the notation in bertsekas.
eps <- .0001                      # tolerance for convergence
k <- 0                            # iteration number
J <- rep(0, length(S))            # initially, start the value iteration with zero cost
J.prev <- rep(eps + 1, length(S)) # to check convergence

while (!all(near(J - J.prev, 0, tol=eps))) {
    k <- k + 1
    J.prev <- J
    for (j in 1:length(S)) {
        J[j] <- FJ(j)
        if (near(k %% 100, 0)) { # print a message to check progress
            message("k = ", k, " J(", j, ") = ", J[j])
        }
    }
}




    
        
