## value iteration for dynamic block size
library(tidyverse)

## state space
X <- seq(0, 200, by=1)  # size of mempool
Y <- seq(0, 25, by=1)   # number of arrivals during prev period

## control is the block size in MB
## 1 MB gets 2000 txs
## 2 MB gets 4000 txs
## and so on
C <- 1:2  # possible controls
K <- 10   # txs per block

alpha <- .9  # discount factor for DP iteration
b <- function(x) 15/(1 + exp(-.05*(x-100))) # tx fee in satoshis per byte
lambda <- .7 # coefficient for AR1 model w_k = lambda*w_{k-1} + xi_k
sigma <- 5   # sd for error term in AR1 model
b0 <- 10     # in the AR1 model

f <- function(x, y, u) {
    xi <- rnorm(1, 0, sigma)
    x.next <- max(0, x - u*K + lambda*y + xi)
    y.next <- max(0, lambda*y + xi)
    list(x.next, y.next)
}

## utility of applying control u in state (x,y)
## really only depends on x_k and not (directly) on y_k
g <- array(data = 0,
           dim = c(length(X), length(C)),
           dimnames = list(statex=as.character(X),
                           control=as.character(C)))
for (x in 1:length(X)) {
    for (u in 1:length(C)) {
        g[x,u] <- b(X[x])*min(C[u]*K, X[x])
    }
}

## applies the DP iteration. note that x,y are the indices
## into the state space
F <- function(x, y) {
    val <- rep(Inf, length(C))   # holds the rev for each control
    for (u in 1:length(C)) {
        cost.to.go <- 0
        for (x.next in X) {
            ## prob of going from X[x],Y[y] to x.next,y.next
            ## see note on computation of cost.to.go
            xi <- x.next - X[x] + min(C[u]*K, X[x]) - lambda*Y[y]
            ## what is the probability that we see this value for xi?
            p <- pnorm(xi+.5, mean=b0, sd=sigma) -
                pnorm(xi-.5, mean=b0, sd=sigma)
            y.next <- min(max(0, lambda*Y[y] + xi), Y[length(Y)])
            cost.to.go <- cost.to.go + p*V[x.next+1,y.next+1]
            ## x.next, y.next are states, but we need to index into V
            val[u] <- g[x,u] + alpha*cost.to.go
        }
    }
    val.star <- max(val)
    ctrl <- which(near(val.star, val))  # could be a tie
    list(val.star, ctrl[1])             # so just choose 1st entry
}
## note on computation of cost.to.go.
## first note that it's really value-to-go. Now, we need to take
## expectation wrt the random component xi. We want to know the
## probability of going from (x,y) to (x.next,y.next).
## from the system equation for x,
## x.next = x - u + lambda*y + xi
## this means that xi has to be equal to
## xi = 
## what is the probability that xi took on this value?
## xi ~ Normal(b0, sigma), so an estimate is provided by
## pnorm(xi+.5, mean=b0, sd=sigma) - pnorm(xi-.5, mean=b0, sd=sigma)
    
## initialization
eps <- .1     # tolerance for convergence
k <- 0          # iteration number
V <- matrix(0, nrow=length(X), ncol=length(Y)) # initial values
V.prev <- matrix(eps + 1, nrow=length(X), ncol=length(Y))
policy <- matrix(NA, nrow=length(X), ncol=length(Y))

## value iteration
while (!all(near(V - V.prev, 0, tol=eps))) {
    k <- k + 1
    V.prev <- V
    for (x in 1:length(X)) {
        for (y in 1:length(Y)) {
            lst <- F(x, y) # pass the indices
            V[x,y] <- lst[[1]]
            policy[x,y] <- lst[[2]]
        }
    }
    if (near(k %% 25, 0)) { # print a message to check progress
        message("k = ", k, "sum(V-V.prev) = ", sum(V-V.prev))
    }
}

## diagnostics and visualization of the optimal value function
## and the optimal policy

SS <- expand.grid(x=X, y=Y, KEEP.OUT.ATTRS = TRUE)
SS$v <- as.vector(V)  # unfolds the matrix column-wise
SS$policy <- as.vector(policy)

ggplot(SS, aes(x=x,y=y)) +
    geom_raster(aes(fill=policy))


library(scatterplot3d)  # an alternative to ggplot2
with(SS, scatterplot3d(x=x, y=y, z=policy, pch=20))
