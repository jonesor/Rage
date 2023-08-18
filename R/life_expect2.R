# Function for mean lifespan. I already wrote essentially the same function in 
# the exactLTRE package, but I want to add the mixing distribution option.
mean_lifespan<- function(Umat, mixdist=NULL, start=NULL){
  # quick check that Umat is square:
  if (dim(Umat)[1]!=dim(Umat)[2]){
    warning('Umat and Fmat are not square matrices.')
  }
  Nclasses<- dim(Umat)[1]
  
  # You cannot use both mixdist and start.
  if (!is.null(mixdist) && !is.null(start)){
    stop('You cannot apply the mixing distribution and also specify a starting state. The mixing distribution defines how you want the function to average over all possible starting states.')
  }
  
  ## Calculate Ex(R | current state)
  N<- exactLTRE::fundamental_matrix(Umat)
  expLCond_z<- rep(1,Nclasses)%*%N
  
  if(!is.null(mixdist)){
    expL<- expLCond_z%*%mixdist
    return(expL)
  } else{
    return(expLCond_z)
  }
  
  if(!is.null(start)){
    return(expLCond_z[start])
  }
}

# Calculate the variance in lifespan:
# note: this calculates the variance in the number of time steps!
var_lifespan<- function(Umat, mixdist=NULL){
  # quick check that Umat is square:
  if (dim(Umat)[1]!=dim(Umat)[2]){
    warning('Umat and Fmat are not square matrices.')
  }
  Nclasses<- dim(Umat)[1]
  
  # calculate the fundamental matrix
  N<- exactLTRE::fundamental_matrix(Umat)
  
  ## Calculate Ex(R | current state)
  expLCond_z<- mean_lifespan(Umat, mixdist = NULL)
  
  ## Var(L | current state) using eqn. 5.12 from Hal's book:
  eT<- matrix(data=1, ncol=Nclasses, nrow=1) # column vector of 1's
  varLCond_z<- eT %*% (2*N%*%N - N) - (expLCond_z)^2
  
  if(is.null(mixdist)){
    return(varLCond_z)
  } else{
    # variance in LRO due to differences along trajectories:
    varL_within<- varLCond_z %*% mixdist 
    # variance in LRO due to differences among starting states:
    varL_between<- t(mixdist)%*%t(expLCond_z^2) - (t(mixdist)%*%t(expLCond_z))^2
    # total variance in lifespan, given the mixing distribution:
    varL<- varL_within + varL_between
    return(varL)
  }
}