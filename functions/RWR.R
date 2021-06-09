RWR <- function(M, p_0, r=0.8 , prop=FALSE, scale = F) {
  
  require("Matrix")

  # use Network propagation when prop=TRUE
  if(prop) {
    w<-Matrix::colSums(M)
    # Degree weighted adjacency matrix
    # define weight matrix, w(i,j) = w(j,i) = sqrt(ki*kj)
    
    W<-matrix(rep(1,nrow(M)*ncol(M)), nrow(M), ncol(M))
    for (i in 1:nrow(M)) {
      for (j in 1:ncol(M)) {
        v<-sqrt(w[i]*w[j])
        W[i,j]<-v
        W[j,i]<-v
      }
    }
    
    # normalise the adj matrix by weight
    W<-M/W
    
    # in the case of non propagation
  } else if(scale) {
    # Convert M to column-normalized adjacency matrix (transition matrix, depend only on the outgoing nodes)
    W=scale(M,center=F,scale=Matrix::colSums(M))
  }
  else{W = M}
  
  # Assign equal probabilities to seed nodes
  
  p_0<-p_0/sum(p_0)
  p_t<-p_0
  # Iterate till convergance is met
  converge = F
  D_t = c(1000, 100) # just some initial large values
  i=1
  while (!converge ) {
    # Calculate new proabalities
    p_tx <- (1-r) * W %*% p_t + r * p_0
    # Check convergance
    D_tx = norm(p_tx-p_t)
    D_t = append(D_t, D_tx)
    print (D_tx)
   
     # convergence is taken if the values don't change, or change very minimally (ratio test)
    if(abs(D_t[i+2]-D_t[i+1])==0){
      converge = T
    } else {
       conv_ratio = log10(abs(D_t[i+2]-D_t[i+1]))/log10(abs(D_t[i+1]-D_t[i]))
    if ( round(conv_ratio, 2) %in% c(0.99,1,1.01) ) {
      converge = T
    } else {
      # converge if lopped for over a hundred times
      if(i >100)
        converge = T
    } #else 
    #f(abs(D_t[i+2]-D_t[i+1]) < 1e-23){
    #  converge = T
    #}
      
       }
    #else{
      #D_t = D_tx
    #}
    p_t<-p_tx
    i= i+1
  }
  
  return(p_tx)
}
