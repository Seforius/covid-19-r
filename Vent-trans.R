library('lpSolve')

#Given T and N we can set up the objective function as below
v <- function(T, N) {
  return(c(replicate((T+1)*N^2,2),replicate((T+1)*N,1e-04), replicate((T+1)*2*N,0)))
  }

#User to set T: Our time periods will start at 0 and end at time T (so T+1). 
#User to set N: the number of nodes, in our case, states of the US. o

T_periods <- 3;
N_nodes <- 7;

#User sets the initial availability of respirators nodes (length = N_nodes)
init_resp <- c(0, 0, 0, 5, 5, 1, 1);

#User sets the required ventilators in each time period at each node
req_vent = matrix(list(), nrow=T_periods+1, ncol=N_nodes) 
req_vent[[1,1]] <- 2;
req_vent[[1,2]] <- 2;
req_vent[[1,3]] <- 1;
req_vent[[1,4]] <- 0;
req_vent[[1,5]] <- 0;
req_vent[[1,6]] <- 2;
req_vent[[1,7]] <- 2;
req_vent[[2,1]] <- 0;
req_vent[[2,2]] <- 2;
req_vent[[2,3]] <- 0;
req_vent[[2,4]] <- 0;
req_vent[[2,5]] <- 0;
req_vent[[2,6]] <- 1;
req_vent[[2,7]] <- 1;
req_vent[[3,1]] <- 0;
req_vent[[3,2]] <- 0;
req_vent[[3,3]] <- 0;
req_vent[[3,4]] <- 4;
req_vent[[3,5]] <- 3;
req_vent[[3,6]] <- 1;
req_vent[[3,7]] <- 1;
req_vent[[4,1]] <- 0;
req_vent[[4,2]] <- 0;
req_vent[[4,3]] <- 0;
req_vent[[4,4]] <- 0;
req_vent[[4,5]] <- 0;
req_vent[[4,6]] <- 0;
req_vent[[4,7]] <- 0;

#Given T and N, we can set up the constraint matrix

zer<-function(x,y)
  return(matrix( rep( 0, len=x*y), nrow = x))

bdiag<-function(e,r) {# e is the block, r the number of repetitions by row
  z = zer(nrow(e), ncol(e)) # zero block of same size as e block
   
 # print(z)
  mr<-matrix(0, ncol = ncol(e)*r, nrow = 0) # initial empty row
  for(i in 1:r) {
    m<-matrix(0, nrow = nrow(e), ncol = 0) # initial empty column
    for(j in 1:r) {
      if (j==i) {
        m <-cbind(m,e)
      #  print(m)
      }
      else {
        m<-cbind(m,z)
      #  print(m)
      }
    }
    mr<-rbind(mr,m)
  }
  return(mr)
}

bsdiag<-function(e,r) {# e is the block, r the number of repetitions by row
  z = zer(nrow(e), ncol(e)) # zero block of same size as e block
 # print(z)
  mr<-matrix(0, ncol = ncol(e)*r, nrow = 0) # initial empty row
  for(i in 1:r) {
    m<-matrix(0, nrow = nrow(e), ncol = 0) # initial empty column
    for(j in 1:r) {
      if (j==i+1) {
        m <-cbind(m,e)
       # print(m)
      }
      else {
        m<-cbind(m,z)
       # print(m)
      }
    }
    mr<-rbind(mr,m)
  }
  return(mr)
}

N = 2
T = 3
e <- matrix(rep(-1,len=N), nrow=1)
print(e)
f <- bdiag(bdiag(e,N),T+1)
#print(f)
h <- bdiag(diag(N), T) + bsdiag(-diag(N),T)
#print(h)
e2<-matrix(0, nrow = N, ncol = 0) # initial empty column
for (i in 1:N) {
  e2 <- cbind(e2, diag(N))
}



#print(e2p)
g <- bdiag(e2, T) + bsdiag(bdiag(e, N), T)
print(g)
print(ncol(g))
print(ncol(h))

#End of constraint matrix set up


C <- c(30, 40, 80)

A <- matrix(c(1, 1, -10,
               4, 3, -20,
               1, 0, -2,
               1, 1, 0), nrow=4, byrow=TRUE)
B <- c(500, 200, 100, 1000)
constranints_direction  <- c("<=", "<=", "<=", ">=")


C <- c(30, 40, 80)

A <- matrix(c(1, 1, -10,
               4, 3, -20,
               1, 0, -2,
               1, 1, 0), nrow=4, byrow=TRUE)
B <- c(500, 200, 100, 1000)
constranints_direction  <- c("<=", "<=", "<=", ">=")
optimum <-  lp(direction="min",
                objective.in = C,
                const.mat = A,
                const.dir = constranints_direction,
                const.rhs = B,
                all.int = T)
print(optimum$status)
best_sol <- optimum$solution
names(best_sol) <- c("x_4p", "x_3p", "x_w")
print(best_sol)
print(paste("Total cost: ", optimum$objval, sep=""))
