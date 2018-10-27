## value iteration for dynamic block size
library(tidyverse)

## state space
X <- seq(0, 100000, by=100)  # size of mempool
Y <- seq(0, 6000, by=100)   # number of arrivals during prev period
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
b <- function(x) 30/(1 + exp(-.0001*(x-50000)))  # tx fee in satoshis per byte
## tx fees were estimated from data available at
## https://www.blockchain.com/charts/mempool-count
## and
## https://bitcoinfees.earn.com/
## a typical tx is about 500 bytes
## a typical tx count in the mempool is about 10,000 txs
## b(10000) returns 11.3 satoshis per byte
## 11.3 * 500 * 1e-8 * 7000 = .4 or about 40 cents for a typical tx fee

lambda <- .7 # coefficient for AR1 model w_k = lambda*w_{k-1} + xi_k
sigma <- 600 # sd for error term in AR1 model
b0 <- 900           # intercept in the AR1 model

## system evolution
## we don't use this function in the value iteration algorithm,
## but it is helpful to understand the underlying process.
f <- function(x, y, u) {
    xi <- rnorm(1, 0, sigma)
    y.next <- max(0, b0 + lambda*y + xi)
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
## x.next = x - u + (b0 + lambda*y + xi)
## and
## y.next = b0 + lambda*y + xi
## During value iteration we will know x,y,x.next,y.next, and u. We
## compute the probability of observing the associated value of xi.
## Note that we can use either system equation to solve for xi.
## Wait! Is that true? The numer of arrivals y is independent of
## the control that is used, but the size of the mempool x is not.
## We want the cost-to-go to depend on the control that is used.
## Furthermore, y and y.next are correlated. It's like this:
## we observe x, y, and y.next and we ask what is the probability
## of going to x.next? xi would have to take on a certain value.
## Let's use the system equation for x.next 
## xi = x.next - x + u - b0 - lambda*y
## what is the probability that xi takes on this value?
## xi ~ Normal(0, sigma), so an estimate is provided by
## pnorm(xi+.5, mean=0, sd=sigma) - pnorm(xi-.5, mean=0, sd=sigma)
P <- array(data = NA,
           dim = c(length(X), length(Y), length(X), length(C)),
           dimnames = list(statex=as.character(X),
                           statey=as.character(Y),
                           statexnext=as.character(X),
                           control=as.character(C)))

library(parallel)

num.cores <- 8
cl <- makeCluster(num.cores, type="FORK")
parSapply(cl, 1:length(X),
          function(x) {

              for (y in 1:length(Y)) {
                  for (x.next in 1:length(X)) {
                      for (u in 1:length(C)) {
                          xi <- X[x.next] - X[x] + C[u]*K - lambda*Y[y]
                          P[x,y,x.next,u] <<- pnorm(xi+.5, mean=0, sd=sigma) -
                              pnorm(xi-.5, mean=0, sd=sigma)
                      }
                  }
              }
              
          })
stopCluster(cl)


require(parallel)
require(doParallel)
num.cores <- 8
cl <- makeCluster(num.cores)
registerDoParallel(num.cores)
P <- foreach(x = 1:length(X), .combine='rbind') %dopar% {
    for (y in 1:length(Y)) {
        for (x.next in 1:length(X)) {
            for (u in 1:length(C)) {
                xi <- X[x.next] - X[x] + C[u]*K - lambda*Y[y]
                P[x,y,x.next,u] <- pnorm(xi+.5, mean=0, sd=sigma) -
                    pnorm(xi-.5, mean=0, sd=sigma)
            }
        }
    }
}
stopCluster(cl)



## initialization
eps <- .01     # tolerance for convergence
k <- 0          # iteration number
V <- matrix(0, nrow=length(X), ncol=length(Y)) # initial values
V.prev <- matrix(eps + 1, nrow=length(X), ncol=length(Y))
policy <- matrix(NA, nrow=length(X), ncol=length(Y))

require(parallel)
require(doParallel)
num.cores <- 28
cl <- makeCluster(num.cores)
registerDoParallel(num.cores)

## value iteration
while (!all(near(V - V.prev, 0, tol=eps))) {
    k <- k + 1
    V.prev <- V
    V <- foreach(x = X, .combine='rbind') %:% 
        foreach(y = Y, .combine='c') %dopar% {
            x.idx <- match(x, X)
            y.idx <- match(y, Y)

            val <- numeric(length(C))   # holds the rev for each control
            for (u in 1:length(C)) {

                cost.to.go <- 0
                ## for each possible x.next, compute the associated y.next
                for (x.next in 1:length(X)) {
                    y.next <- X[x.next] - x + C[u]*K  # y.next is a state not an index
                    y.next <- ifelse(y.next<0, 0, y.next)
                    y.next <- ifelse(y.next>max(Y), max(Y), y.next)
                    y.next.idx <- match(y.next, Y)       # locate the index
                    
                    cost.to.go <- cost.to.go + P[x.idx,x.next,y.idx,u]*V[x.next,y.next.idx]
                }
                val[u] <- g[x.idx,u] + cost.to.go
            }
            V[x.idx,y.idx] <- max(val)
            ctrl <- which(near(max(val), val))[1]  # could be a tie
            policy[x.idx,y.idx] <- C[ctrl]
        }
    message("k = ", k, "sum(V-V.prev) = ", sum(V-V.prev))
}

stopCluster(cl)

## save results to disk
save.image("dp-blocksize.RData")

if (!interactive()) stop()

## diagnostics and visualization of the optimal value function
## and the optimal policy
SS <- expand.grid(x=X, y=Y, KEEP.OUT.ATTRS = TRUE)
SS$v <- as.vector(V)  # unfolds the matrix column-wise
SS$policy <- as.vector(policy)

ggplot(SS, aes(x=x,y=y)) + geom_raster(aes(fill=policy))
