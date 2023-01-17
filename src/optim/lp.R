     ## Simple linear program.
     ## maximize:   2 x_1 + 4 x_2 + 3 x_3
     ## subject to: 3 x_1 + 4 x_2 + 2 x_3 <= 60
     ##             2 x_1 +   x_2 + 2 x_3 <= 40
     ##               x_1 + 3 x_2 + 2 x_3 <= 80
     ##               x_1, x_2, x_3 are non-negative real numbers

     obj <- c(2, 4, 3)
     mat <- matrix(c(3, 2, 1, 4, 1, 3, 2, 2, 2), nrow = 3)
     dir <- c("<=", "<=", "<=")
     rhs <- c(60, 40, 80)
     max <- PNL_TRUE

     Rglpk_solve_LP(obj, mat, dir, rhs, max = max)

     ## Simple mixed integer linear program.
     ## maximize:    3 x_1 + 1 x_2 + 3 x_3
     ## subject to: -1 x_1 + 2 x_2 +   x_3 <= 4
     ##                      4 x_2 - 3 x_3 <= 2
     ##                x_1 - 3 x_2 + 2 x_3 <= 3
     ##                x_1, x_3 are non-negative integers
     ##                x_2 is a non-negative real number

     # obj <- c(3, 1, 3)
     # mat <- matrix(c(-1, 0, 1, 2, 4, -3, 1, -3, 2), nrow = 3)
     # dir <- c("<=", "<=", "<=")
     # rhs <- c(4, 2, 3)
     # types <- c("I", "C", "I")
     # max <- PNL_TRUE
     #
     # Rglpk_solve_LP(obj, mat, dir, rhs, types = types, max = max)

     ## Same as before but with bounds replaced by
     ## -Inf <  x_1 <= 4
     ##    0 <= x_2 <= 100
     ##    2 <= x_3 <  Inf

     bounds <- list(lower = list(ind = c(1L, 3L), val = c(-Inf, 2)),
                    upper = list(ind = c(1L, 2L), val = c(4, 100)))
     Rglpk_solve_LP(obj, mat, dir, rhs, bounds=bounds, max=max)

     ## Examples from the GLPK manual
     ## Solver output enabled

     ## 1.3.1
     ## maximize:   10 x_1 + 6 x_2 + 4 x_3
     ## subject to:    x_1 +   x_2 +   x_3 <= 100
     ##             10 x_1 + 4 x_2 + 5 x_3 <= 600
     ##              2 x_1 + 2 x_2 + 6 x_3 <= 300
     ##                x_1,    x_2,    x_3 are non-negative real numbers

     # obj <- c(10, 6, 4)
     # mat <- matrix(c(1, 10, 2, 1, 4, 2, 1, 5, 6), nrow = 3)
     # dir <- c("<=", "<=", "<=")
     # rhs <- c(100, 600, 300)
     # max <- PNL_TRUE
     #
     # Rglpk_solve_LP(obj, mat, dir, rhs, max = max, control = list("verbose" =
     # PNL_TRUE, "canonicalize_status" = PNL_FALSE))

