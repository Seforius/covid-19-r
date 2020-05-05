library('lpSolve')

#Given T and N we can set up the objective function as below
v <- function(T, N) {
  return(c(replicate((T+1)*N^2,2),replicate((T+1)*N,1e-04), replicate((T+1)*2*N,0)))
  }

#function to build a zero matrix that is x by y
zer<-function(x,y)
  return(matrix( rep( 0, len=x*y), nrow = x, ncol = y))

#function to build a block diagonal matrix. e is the block, r the repetitions.

bdiag<-function(e,r) {# e is the block, r the number of repetitions by row
  z = zer(nrow(e), ncol(e)) # zero block of same size as e block
  print(z)
  mr<-matrix(0, ncol = ncol(e)*r, nrow = 0) # initial empty row
  for(i in 1:r) {
    m<-matrix(0, nrow = nrow(e), ncol = 0) # initial empty column
    for(j in 1:r) {
      if (j==i) {
        m <-cbind(m,e)
        print(m)
      }
      else {
        m<-cbind(m,z)
        print(m)
      }
    }
    mr<-rbind(mr,m)
  }
  return(mr)
}

#Function to build a matrix that is almost block diagonal, shiofted to the right by one block.
bsdiag<-function(e,r) {# e is the block, r the number of repetitions by row
  z = zer(nrow(e), ncol(e)) # zero block of same size as e block
  print(z)
  mr<-matrix(0, ncol = ncol(e)*r, nrow = 0) # initial empty row
  for(i in 1:r) {
    m<-matrix(0, nrow = nrow(e), ncol = 0) # initial empty column
    for(j in 1:r) {
      if (j==i+1) {
        m <-cbind(m,e)
        print(m)
      }
      else {
        m<-cbind(m,z)
        print(m)
      }
    }
    mr<-rbind(mr,m)
  }
  return(mr)
}

#mat function below will build the constraint matrix using the functions above that develop bd matrics.

mat <-function(T, N) {
  a <- zer((T+1)*N, (T+1)*N^2)
  b <- diag((T+1)*N)
  c <- zer((T+1)*N, (T+1)*N)
  d <- diag(N)
  e <- matrix(rep(-1,len=N), nrow=1)
  e2<-matrix(0, nrow = N, ncol = 0) # initial empty column
  for (i in 1:N) {
    e2 <- cbind(e2, diag(N))
  }
  

  
  f <- bdiag(bdiag(e,N),T+1) 
  g <- bdiag(e2, T) + bsdiag(bdiag(e, N), T)
  g2 <- rbind(zer(N*(T-1),N^2), bdiag(e, N))
  g <-cbind(g,g2)
  h <- bdiag(diag(N), T) + bsdiag(-diag(N),T)
  h2 <- rbind(zer(N*(T-1),N),-diag(N))
  h <-cbind(h,h2)
  i <- zer(N,(T+1)*N^2)
  j <- zer(N,(T+1)*N)
  k <- cbind(diag(N), zer(N,T*N))
  l <- zer(N*T,(T+1)*N)
  
  
  r1 <- cbind(a,b,b,c)
  print(ncol(r1))
  r2 <- cbind(a,c,-b,b)
  print(ncol(r2))
  r3 <- cbind(f,c,c,b)
  print(ncol(r3))
  r4 <- cbind(g,l,h,l)
  print(ncol(r4))
  r5 <- cbind(i,j,k,j)
  print(ncol(r5))
  
  return(rbind(r1,r2,r3,r4,r5))
}

#User to set T: Our time periods will start at 0 and end at time T (so T+1). 
#User to set N: the number of nodes, in our case, states of the US. o

T_periods <- 3;
N_nodes <- 7;

#User sets the initial availability of respirators nodes (length = N_nodes)
init_resp <- c(0, 0, 0, 5, 5, 1, 1);

#User sets the required ventilators in each time period at each node (list by fixing t first)
req_vent <- c(2,2,1,0,0,2,2,0,2, 0,0,0,1,1,0,0,0,4,3,1,1,0,0,0,0,0,0,0) 

#Given T and N, we can set up the objective coefficients 

C <- v(T_periods, N_nodes) 

#Then set up the constrain matrix: 

A <- mat(T_periods, N_nodes)

#Then fill up the right hand side
B <- c(req_vent, req_vent, rep(0, ((2*T_periods+1)*N_nodes)), init_resp) 

#Lastly fill in the constraint directions
constraints_direction  <- c("=")
for(i in 2:((2*(T_periods+1))*N_nodes))
{
  constraints_direction <- c(constraints_direction, "=")
}

for(i in 1:(((T_periods+1))*N_nodes))
{
  constraints_direction <- c(constraints_direction, ">=")
}

for(i in 1:(((T_periods+1))*N_nodes))
{
  constraints_direction <- c(constraints_direction, "=")
}

optimum <-  lp(direction="min",
                objective.in = C,
                const.mat = A,
                const.dir = constraints_direction,
                const.rhs = B,
                all.int = F)
print(optimum$status)
best_sol <- optimum$solution
print(best_sol)
print(paste("Total cost: ", optimum$objval, sep=""))
