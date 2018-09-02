## value iteration for dynamic block size
library(tidyverse)

## state space
X <- seq(0, 100000, by=10)  # size of mempool
Y <- seq(0, 6000, by=10)   # number of arrivals during prev period
## these numbers were arrived at by looking at the history
## from aug '17 to aug '18 (which includes the dec/jan busy period)
## of the
## mempool tx count (num of txs waiting to be confirmed), and the
## txs arrival rate (num of txs added to the mempool per second)
## source: https://www.blockchain.com/charts

## control is the block size in MB
## 1 MB gets 2000 txs
## 2 MB gets 4000 txs
## and so on
C <- 1:3   # possible controls
K <- 2000  # txs per unit of control

alpha <- .9          # discount factor for DP iteration
b <- function(x) 30/(1 + exp(-.0001*(x-15000)))  # tx fee in satoshis per byte
## tx fees were estimated from data available at
## https://www.blockchain.com/charts/mempool-count
## and
## https://bitcoinfees.earn.com/
## a typical tx is about 500 bytes
## a typical tx count in the mempool is about 10,000 txs
## b(10,000) returns 11.3 satoshis per byte
## 11.3 * 500 * 1e-8 * 7000 = .4 or about 40 cents for a typical tx fee

lambda <- .7         # coefficient for AR1 model w_k = lambda*w_{k-1} + xi_k
sigma <- sqrt(4.7e6) # sd for error term in AR1 model
b0 <- 0              # intercept in the AR1 model

## system evolution
## we don't use this function in the value iteration algorithm,
## but it is helpful to understand the underlying process.
f <- function(x, y, u) {
    xi <- rnorm(1, b0, sigma)
    y.next <- max(0, lambda*y + xi)
    x.next <- max(0, x - u*K + y.next)
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
## what is the probability that xi takes on this value?
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
            val <- rep(Inf, length(C))   # holds the rev for each control
            for (u in 1:length(C)) {
                ## possible values of x.next for each possible y.next
                x.next <- X[x] - C[u]*K + Y
                x.next <- ifelse(x.next<0, 0, x.next)
                x.next <- ifelse(x.next>max(X), max(X), x.next)
                ## note that x.next is a vector of states, not indices
                ## so we need to determine the associated indices
                x.next.idx <- match(x.next, X)
                val[u] <- g[x,u] +
                    alpha * sum(P[y,1:length(Y)] * diag(V[x.next.idx,1:length(Y)]))
            }
            V[x,y] <- max(val)
            ctrl <- which(near(max(val), val))[1]  # could be a tie
            policy[x,y] <- C[ctrl]
        }
    }
    if (near(k %% 10, 0)) { # print a message to check progress
        message("k = ", k, "sum(V-V.prev) = ", sum(V-V.prev))
    }
}

## save results to disk
save.image("dp-blocksize.RData")
 
## diagnostics and visualization of the optimal value function
## and the optimal policy
SS <- expand.grid(x=X, y=Y, KEEP.OUT.ATTRS = TRUE)
SS$v <- as.vector(V)  # unfolds the matrix column-wise
SS$policy <- as.vector(policy)

ggplot(SS, aes(x=x,y=y)) + geom_raster(aes(fill=policy))

