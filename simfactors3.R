### Code ported from SAS IML program in Coughlin et al. 2013 (SESUG 2013)
###
### Scalars follow the SAS/IML code for now. Will eventually modify
### some to be more flexible, including pulling the necessary data out
### of a factor loading matrix.

### Section includes other preliminary stuff, too.

## Set random seed to something
seed <- 123

set.seed(seed)

## Matrix to enable selection of communalities to be done without if-then statements
comm_select <- cbind(c(.2999999, .6999999, .2999999), c(.15, .15, .55))

## Scalars
nvars <- 10 # Number of variables (rows)

nfactors <- 3 # Number of factors (columns)

d_frac <- 0 # Proportion of variables that are dichotomous; uses
## the product of this and nvars is rounded to an integer

commun_type <- 3 # Communalities 1 = low; 2 = wide; 3 = high

## For commun_type modify at some point so that it can extract
## communalities from a loading matrix

replicat <- 10 # Number of samples drawn from each population

numobs <- 150 # Number of observations per sample drawn

### Create various matrices, etc. In original SAS/IML code there are
### several steps used to generate. With the exception of what looks
### like the need to set up matrices in advance in SAS/IML, I will
### preserve most of them for the time being as this makes it easier
### to debug the code

##samples <- vector("list", 10)
##for (ii in 1:10) {

## Generate high, low, or wide communalities
bp2 <- matrix(runif(nvars ^ 2), ncol = nvars)

bp3 <- (bp2 * comm_select[commun_type, 1]) + comm_select[commun_type, 2]

b1square <- (round(diag(bp3), 1)) * diag(nvars)

B1 <- b1square ^ .5

b3square <- diag(nvars) - b1square

B3 <- b3square ^ .5

### Generate conceptual factor loadings matrix (A1~) based on method
### of Tucker, Koopman, and Linn (1969). There is probably a more
### efficient way to do this, but the method below follows the
### procedure to the letter.
A1tilde <- (matrix(0, nrow = nvars, ncol = nfactors))

for (i in 1:nvars) {

    ## Creates vector that indicates what conceptual loading is in
    ## column 1, 2 ... nfactors for loop below, which is where the value is
    ## placed.
    selector <- sample(nfactors, replace = F)

    ## Loop that picks random values for p-1 of the conceptual loadings.
    for (j in 1:nfactors - 1) {

        A1tilde[i, selector[j]] <- sample((0:(nfactors - 1 -
                                              sum(A1tilde[i, ]))), 1, replace = T)
    }

    ## The last conceptual loading should be equal to the highest
    ## possible value (p-1) minus the sum of the existing conceptual
    ## loadings.
    A1tilde[i, selector[nfactors]] <- nfactors - 1 - sum(A1tilde[i, ])
}

x <- matrix(rnorm(A1tilde, mean = 0, sd = 1), nrow = nvars,
            ncol = nfactors)

x2 <- x ^ 2

## Need to check next line in SAS (checked and okay)
d <- matrix(rowSums(x2) ^ -.5, nrow = nvars, ncol = nfactors)

## Need to check next line with SAS (checked and okay)
cvec <- matrix(round((runif(nfactors) * .2999999) + .65, 1), nrow = 1,
               ncol = nfactors)

c <- matrix(1, nrow = nvars, ncol = 1) %*% cvec

c2 <- c ^ 2

ones <- matrix(1, nrow = nvars, ncol = nfactors)

y <- A1tilde * c + d * x * ((ones - c2) ^ .5)

k <- .2

## Tricky bit of code to create the z-matrix, namely as I'm doing it via vectorization and not a loop!
## SAS code gives same results, so all is okay!
z <- matrix(((1 + k) *
             (matrix(y, ncol = 1)) *
             ((matrix(y, ncol = 1)) + abs(matrix(y, ncol = 1)) + k)) / ((2 + k) *
                                                                        (abs(matrix(y, ncol = 1)) + k)), ncol = nfactors)

## Final bit where it's all put together to create the population correlation matrix R
z2 <- z ^ 2

g <- matrix(rowSums(z2) ^ -.5, nrow = nvars, ncol = nfactors)

A1star <- g * z

A1 <- B1 %*% A1star

A3star <- diag(1, nvars)

A3 <- B3 * A3star

R <- A1 %*% t(A1) + A3 %*% t(A3)

##samples[ii] <- list(R)


##}

### Cholesky root approach to generating random data
### Code based on https://www.r-bloggers.com/simulating-random-multivariate-correlated-data-continuous-variables/
samples <- vector("list", replicat)

i <- 1

while(i <= replicat) {

    U <- t(chol(R))

    X <- U %*% (matrix(rnorm(nvars * numobs, 0, 1), nrow = nvars, ncol = numobs))

    newX <- t(X)

    if(!is.na(d_frac) & d_frac > 0)

        { newX[, 1:(round(nvars * d_frac))] <- as.numeric(newX[, 1:(round(nvars * d_frac))] > 0) }

    samples[i] <- list(newX)

    i <- i + 1
}

#### Code below seems to be generating random values from correlation matrix

## XX <- matrix(rnorm(nvars * numobs, 0, 1), nrow = numobs, ncol = nvars)

## UU <- eigen(R)$vectors
## DD <- diag(eigen(R)$values)
## FF <- UU%*%DD^.5
## ZZZ<-XX %*% t(FF)

fax <- function(x) {fa(x,nfactor=4)}
