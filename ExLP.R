library('lpSolve')

C <- c(2, 5)

A <- matrix(c(3, 4,
               9, 7,
               1, 0,
               0, 1), nrow=4, byrow=TRUE)
B <- c(200, 300, 20, 0)
constraints_direction  <- c("<=", "<=", ">=", ">=")
optimum <-  lp(direction="max",
                objective.in = C,
                const.mat = A,
                const.dir = constraints_direction,
                const.rhs = B,
                all.int = T)
print(optimum$status)
best_sol <- optimum$solution
names(best_sol) <- c("x1", "x2")
print(best_sol)
print(paste("Total units: ", optimum$objval, sep=""))

