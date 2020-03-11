I.rule <- "min"
II.fobj <- c(6,8,10,12,2,2,2,1000,1000,1000)

III.Acon <- matrix(c(1,0,0,0,-1,0,0,0,0,0,
                     0,1,0,0,1,-1,0,0,0,0,
                     0,0,1,0,0,1,-1,0,0,0,
                     0,0,0,1,0,0,1,0,0,0,
                     1,0,0,0,0,0,0,0,0,0,
                     0,1,0,0,0,0,0,0,0,0,
                     0,0,1,0,0,0,0,0,0,0,
                     0,0,0,1,0,0,0,0,0,0), 
                   nrow = 8, 
                   byrow = TRUE)

III.Acon
dim(III.Acon)

IV.dir <- c("=",
            "=",
            "=",
            "=",
            "<=",
            "<=",
            "<=",
            "<=")
IV.dir
length(IV.dir)

V.bound <- c(100,
             200,
             150,
             400,
             400,
             400,
             300,
             300)
V.bound
length(V.bound)

library("lpSolve")
case2p2.sol <- lp(I.rule,
                  II.fobj,
                  III.Acon,
                  IV.dir,
                  V.bound)

case2p2.sol$solution

case2p2.sol$objval

