I.rule <- "max"
II.fobj <- c(1,9,1)
III.Acon <- matrix(c(1,2,3,3,2,2), nrow = 2, byrow = TRUE)
IV.dir <- c("<=","<=")
V.bound <- c(9,15)

library("lpSolve")
case1.sol <- lp(I.rule,
                II.fobj,
                III.Acon,
                IV.dir,
                V.bound)

case1.sol$objval

case1.sol$solution
