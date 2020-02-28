I.rule <- "min"
II.fobj <- c(12,12,12,12,2,2,2,2)
III.Acon <- matrix(c(1,0,0,0,-1,0,0,0,
                     0,1,0,0,1,-1,0,0,
                     0,0,1,0,0,1,-1,0,
                     0,0,0,1,0,0,1,-1,
                     1,0,0,0,0,0,0,0,
                     0,1,0,0,0,0,0,0,
                     0,0,1,0,0,0,0,0,
                     0,0,0,1,0,0,0,0), 
                   nrow = 8, 
                   byrow = TRUE)

IV.dir <- c("=",
            "=",
            "=",
            "=",
            "<=",
            "<=",
            "<=",
            "<=")

V.bound <- c(100,
             200,
             150,
             400,
             400,
             400,
             300,
             300)

library("lpSolve")
case2.sol <- lp(I.rule,
                II.fobj,
                III.Acon,
                IV.dir,
                V.bound)

case2.sol$objval

case2.sol$solution
