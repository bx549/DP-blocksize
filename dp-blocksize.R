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
K <- 10   # txs per unit of control

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
blk.size <- sapply(X, function(x) pmin(C*K, x))  # num txs in block
g <- b(X) * t(blk.size)  # element-wise

## transition probabilities.
## When computing the value-to-go, we need to take
## expectation wrt the random component xi. We want to know the
## probability of going from (x,y) to (x.next,y.next) when control u
## is applied. From the system equations for x.next, y.next
## x.next = x - u + (lambda*y + xi)
## and
## y.next = lambda*y + xi
## During value iteration we will know x,y,x.next,y.next, and u. We
## compute the probability of observing the associated value of xi.
## Note that we can use either system equation to solve for xi.
## Let's use the system equation for y.next since the state space
## for y is smaller than for x. Given y and y.next,
## xi = y.next - lambda*y
## what is the probability that xi took on this value?
## xi ~ Normal(b0, sigma), so an estimate is provided by
## pnorm(xi+.5, mean=b0, sd=sigma) - pnorm(xi-.5, mean=b0, sd=sigma)
## Note that the distribution of the random xi does not depend on the
## control u.
P <- array(data = NA,
           dim = c(length(Y), length(Y)),
           dimnames = list(statey=as.character(Y),
                           stateynext=as.character(Y)))
for (y in 1:length(Y)) {
    xi <- Y - lambda*Y[y]
    P[y,] <- pnorm(xi+.5, mean=b0, sd=sigma) -
        pnorm(xi-.5, mean=b0, sd=sigma)
}

## applies the DP iteration.
## note that x,y are indices into the state space, while x.next, y.next
## are state space values.
F <- function(x, y) {
    val <- rep(Inf, length(C))   # holds the rev for each control
    for (u in 1:length(C)) {
        val.to.go <- 0
        for (y.next in Y) {
            ## expectation is wrt xi.
            ## see note on computation of transition probabilities.
            x.next <- min(max(0, X[x] - C[u]*K + y.next), X[length(X)])
            val.to.go <- val.to.go + P[y,y.next+1]*V[x.next+1,y.next+1]
            val[u] <- g[x,u] + alpha*val.to.go
        }
    }
    val.star <- max(val)
    ctrl <- which(near(val.star, val))  # could be a tie
    list(val.star, ctrl[1])             # so just choose 1st entry
}
    
## initialization
eps <- .01     # tolerance for convergence
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

ggplot(SS, aes(x=x,y=y)) + geom_raster(aes(fill=policy))

