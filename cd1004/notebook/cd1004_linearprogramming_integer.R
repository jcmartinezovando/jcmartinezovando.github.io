# Resolucion de problemas con valores enteros

if(!require("lpSolve")){install.packages("lpSolve")}

library("lpSolve")

# Pensemos en el siguiente sistema de programacion lineal:
#
# MAXIMIZAR
# x1 + 9 x2 + x3 
#
# Sujeto a:
# x1 + 2 x2 + 3 x3 <= 9
# 3 x1 + 2 x2 + 2 x3 <= 15

a.rule <- "max"
b.obj <- c(1, 9, 1)
c.con <- matrix (c(1, 2, 3, 3, 2, 2), nrow=2, byrow=TRUE)
d.dir <- c("<=", "<=")
e.rhs <- c(9, 15)

lp.solution = lp(a.rule, b.obj, c.con, d.dir, e.rhs)

lp.solution$solution
lp.solution$objval


# Dense constraint:
#
c.con.d <- matrix(c(rep (1:2,each=3), rep (1:3, 2), t(c.con)), ncol=3)
c.con
c.con.d

lp.solution.d = lp(a.rule, b.obj, c.con.d, d.dir, e.rhs)

lp.solution.d$solution
lp.solution.d$objval

# Sensitiveness
lp.solution.s = lp(a.rule, b.obj, c.con, d.dir, e.rhs, compute.sens=TRUE)
lp.solution.s$sens.coef.from
lp.solution.s$sens.coef.to

lp.solution.s$duals
lp.solution.s$duals.from
lp.solution.s$duals.to


# Integer
lp.solution.int = lp (a.rule, b.obj, c.con, d.dir, e.rhs, int.vec=2)

lp.solution.int$solution
lp.solution.int$objval