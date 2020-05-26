# Solo primer ITEM

library("lpSolve")

a.rule <- "min"

#            Q1  Q2  Q3  Q4
b.fobj <- c( 12, 14, 16, 18)

#                  Q1 Q2 Q3 Q4
c.Acon <- matrix(c(1, 0, 0, 0,     # Mes 1 Demanda:  (1)Q1 + (0)Q2 + (0)Q3 + (0)Q4 = 200
                   0, 1, 0, 0,
                   0, 0, 1, 0,
                   0, 0, 0, 1,
                   1, 0, 0, 0,     # Mes 1 Workinghours:  (1)Q1 + (0)Q2 + (0)Q3 + (0)Q4 <= 300
                   0, 1, 0, 0,
                   0, 0, 1, 0,
                   0, 0, 0, 1), 
                   nrow = 8, 
                   byrow = TRUE)
  
d.dir <- c( "=",      # Mes 1 Demanda:  (1)Q1 + (0)Q2 + (0)Q3 + (0)Q4 = 200
            "=",
            "=",
            "=", 
            "<=",
            "<=",
            "<=",
            "<=")

e.limits <- c( 100,      # Mes 1 Demanda:  (1)Q1 + (0)Q2 + (0)Q3 + (0)Q4 = 200
               200,
               150,
               400,
               400,     # Mes 1 Workinghours:  (1)Q1 + (0)Q2 + (0)Q3 + (0)Q4 <= 300
               400,
               300,
               300)

library("lpSolve")
esta.es.mi.sol <- lp(a.rule,
                     b.fobj,
                     c.Acon,
                     d.dir,
                     e.limits)

# Estrategia optima
esta.es.mi.sol$solution

# Valor de la funcion en la estrategia optima
esta.es.mi.sol$objval
