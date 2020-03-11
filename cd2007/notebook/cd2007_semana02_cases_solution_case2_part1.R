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

rule <- "min"

IV.dir <- c("=",
            "=",
            "=",
            "=",
            "<=",
            "<=",
            "<=",
            "<=")

V <- c(100,
             200,
             150,
             400,
             400,
             400,
             300,
             300)

fobj <- c(12,12,12,12,2,2,2,2)



library("lpSolve")
esta.es.mi.sol <- lp(rule,
                fobj,
                III.Acon,
                IV.dir,
                V.bound)

case2.sol$objval

case2.sol$solution
