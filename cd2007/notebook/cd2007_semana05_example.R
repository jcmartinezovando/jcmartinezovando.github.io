f.obj <- c(1, 9, 1)
f.con <- matrix (c(1, 2, 3, 3, 2, 2), nrow=2, byrow=TRUE)
f.dir <- c("<=", "<=")
f.rhs <- c(9, 15)


library("lpSolve")
example.sol <- lp ("max", f.obj, f.con, f.dir, f.rhs)

example.sol$solution

# Sensitiveness

lp ("max", f.obj, f.con, f.dir, f.rhs, compute.sens=TRUE)$sens.coef.from

lp ("max", f.obj, f.con, f.dir, f.rhs, compute.sens=TRUE)$sens.coef.to

# Duals

lp ("max", f.obj, f.con, f.dir, f.rhs, compute.sens=TRUE)$duals

# Mixed-type

mixtype.sol <- lp ("max", f.obj, f.con, f.dir, f.rhs, int.vec=1:3)

mixtype.sol$solution

# 8-queens problem,
chess.obj <- rep (1, 64)
chess.obj

q8 <- make.q8 ()
q8

chess.dir <- rep (c("=", "<"), c(16, 26))
chess.dir

chess.rhs <- rep (1, 42)
chess.rhs

queens8.sol <- lp ('max',
                    chess.obj, , 
                    chess.dir, 
                    chess.rhs, 
                    dense.const = q8,
                    all.bin=TRUE, 
                    num.bin.solns=3)

queens8.sol$solution