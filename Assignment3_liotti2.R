myBW = function(x, A, B, w, n.iter = 100){
  # Input: 
  # x: T-by-1 observation sequence
  # A: initial estimate for mz-by-mz transition matrix
  # B: initial estimate for mz-by-mx emission matrix
  # w: initial estimate for mz-by-1 initial distribution over Z_1
  # Output MLE of A and B; we do not update w
  # list(A = A, B=B, w = w)
  
  for(i in 1:n.iter){
    update.para = BW.onestep(x, A, B, w)
    A = update.para$A
    B = update.para$B
  }
  return(list(A = A, B = B, w = w))
}

BW.onestep = function(x, A, B, w){
  # Input: 
  # x: T-by-1 observation sequence
  # A: current estimate for mz-by-mz transition matrix
  # B: current estimate for mz-by-mx emission matrix
  # w: current estimate for mz-by-1 initial distribution over Z_1
  # Output the updated parameters
  # para = list(A = A1, B = B1)
  
  # We DO NOT update the initial distribution w
  
  T = length(x)
  mz = nrow(A)
  alp = forward.prob(x, A, B, w)
  beta = backward.prob(x, A, B, w)
  myGamma = array(0, dim=c(mz, mz, T-1))
  
  ###
  ## YOUR CODE: 
  ## Compute gamma_t(i,j), which are stored in myGamma
  ##
  
  for (t in 1:dim(myGamma)[3]) {
    for (i in 1:mz) {
      for (j in 1:mz) {
        myGamma[i,j,t] = alp[t, i] * A[i,j] * B[j, x[t+1]] * beta[t+1, j]
      } 
    }
  }
  
  A = rowSums(myGamma, dims = 2)
  A = A/rowSums(A)
  
  tmp = apply(myGamma, c(1, 3), sum)  # mz-by-(T-1)
  tmp = cbind(tmp, colSums(myGamma[, , T-1]))
  for(l in 1:mx){
    B[, l] = rowSums(tmp[, which(x==l)])
  }
  B = B/rowSums(B)
  return(list(A = A, B = B))
}

forward.prob = function(x, A, B, w){
  
  # Output the forward probability matrix alp 
  # alp: T by mz, (t, i) entry = P(x_{1:t}, Z_t = i)
  
  T = length(x)
  mz = nrow(A)
  alp = matrix(0, T, mz)
  
  # fill in the first row of alp
  alp[1, ] = w*B[, x[1]]
  
  # Recursively compute the remaining rows of alp
  for(t in 2:T){
    tmp = alp[t-1, ] %*% A
    alp[t, ] = tmp * B[, x[t]]
  }
  return(alp)
}

backward.prob = function(x, A, B, w){
  # Output the backward probability matrix beta
  # beta: T by mz, (t, i) entry = P(x_{(t+1):n} | Z_t = i)
  # for t=1, ..., n-1
  
  T = length(x)
  mz = nrow(A)
  beta = matrix(1, T, mz)
  
  # The last row of beta is all 1.
  # Recursively compute the previous rows of beta
  for(t in (T-1):1){
    tmp = as.matrix(beta[t+1, ] * B[, x[t+1]])  # make tmp a column vector
    beta[t, ] = t(A %*% tmp)
  }
  return(beta)
}

myViterbi = function(data, A, B, iniw) {
  delta = matrix(0, nrow = length(data), ncol = 2)
  n = length(data)
  
  for (k in 1:2) {
    delta[1, k] = log(iniw[k]) + log(B[k, data[1]])
  }
  
  for (t in 2:nrow(delta)) {
    for (i in 1:2) {
      j = which.max(delta[t-1, ])
      delta[t, i] = max(delta[(t-1), j] + log(A[j, i])) + log(B[i, data[t]])
    }
  }
  
  myout.z = rep(0, n)
  
  myout.z[n] = which.max(delta[n, ])
  l = 1:2
  for (i in (n-1):1) {
      myout.z[i] = which.max(delta[i, l] + log(A[l, myout.z[i+1]]))
  }
  
  for (i in 1:n) {
    
    if (myout.z[i] == 1) {
      myout.z[i] = "A"
    }
    
    if (myout.z[i] == 2) {
      myout.z[i] = "B"
    }
    
  }
  
  return(myout.z)
} 

data = read.csv("Coding3_HMM_Data.csv")

mz=2; mx=3
ini.A = matrix(1, mz, mz)
ini.A = ini.A/rowSums(ini.A)
ini.B = matrix(1:6, mz, mx)
ini.B = ini.B/rowSums(ini.B)
ini.w = c(1/2, 1/2)

myout = myBW(data$X, ini.A, ini.B, ini.w, n.iter = 100)
myout.Z = myViterbi(data$X, myout$A, myout$B, ini.w)
write.table(myout.Z, file = "Coding3_HMM_Viterbi_Output.txt", 
            row.names = FALSE, col.names = FALSE)

########Checking accuracy
library(HMM)
hmm0 = initHMM(c("A", "B"), c(1, 2, 3), 
              startProbs = ini.w,
              transProbs = ini.A, emissionProbs = ini.B)

true.out = baumWelch(hmm0, data$X, maxIterations=100, pseudoCount=0)
true.out$hmm

true.viterbi = viterbi(true.out$hmm, data$X)
#write.table(true.viterbi, file = "Coding3_HMM_True_Viterbi_Output.txt", 
#            row.names = FALSE, col.names = FALSE)
sum(true.viterbi != myout.Z)
