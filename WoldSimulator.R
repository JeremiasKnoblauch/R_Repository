#THIS IS CHANGED


#Generate stationary processes as you wish, and obtain the corresponding TRUE Parameters
#
#   ERROR type number code for k = 1,2,3, ...
#   1   normal(0,1)
#   k   t-distr df=k; k>2
#
#   COEFFICIENT sequence code l \in (0,1]
#   mode=.. c_j =...
#   0       c_j = 1/(j^l) for l<1, mode = 0
#   1       c_j = log_2(1+j^l) for l<1, mode = 1
#   2       c_j = a^j, mode=2
#   randomPeaksBoolean=TRUE:  
#           generate random Peaks whose size
#           is uniformly distr. on (-randomPeaksRange, randomPeaksrange)
#           and whose position is randomly designed, also uniformly over 1:n.
#           The frequency of thse random peaks is determined by randomPeaksFrequency
#
#   FUNCTIONS and what they do
#   errorGenerator(n,T,k)                 ->  using errors of code k, errors are produced for
#                                             a process with T observations and n coefficients on the errors (so n+T 
#                                             errors are produced)  
#   coefficientGenerator(n,l, mode, randomPeaksBoolean, randomPeaksFrequency=NULL, randomPeaksRange=NULL, a=NULL) 
#                                         ->  using coefficients corresponding to the code, a sequence of length n
#                                             is produced of these coefficients. See 'mode' and 'randomPeaksBoolean' above
#   WoldGenerator(n,T,k,coefSequence)     ->  using n,T,k, errorGenerator is called. The readily produced
#                                             sequence of coefficients (coefsequence) is then used to generate
#   trueGammas(coefSequence,p,k)          ->  calculates the covariance matrix of the asymptotic normal (Gamma_p),
#                                             the first p autocovariances (gamma_p) such that alpha<-Gamma_p^{-1}%*%gamma_p,
#                                             the process/long run variance (gamma_0) and the error term variance (sigma2).
#                                             All are stored in a list and accessible with these names.
#   trueAlphas(trueGammaList)             ->  Takes the output of trueGammas(.,.,.) and processes it into the true alpha-projection
#                                             coefficients via Hamilton's formula: alpha<-Gamma_p^{-1}%*%gamma_p.
#                                         
#NOTE 1: I do not need a burn-in because I work directly with the error terms (rather than the x_t)
#NOTE 2: This is different with AR-generators (which I might introduce later)
#
#___________________________________________________________________________________________________________________

#clear objects
#rm(list=ls())
#gc()

#generator for errors, return as E (normals or Student-t)
errorGenerator<-function(n,T,k)
{
  if(k==1) #use standard normal
  {
    E<-matrix(rnorm(n+T), ncol=1)
  }
  else #use t with k degrees of freedom
  {
    E<-matrix(rt(n+T, k), ncol=1)
  }
  return(E)
}

#generator for coefficients, return as C
coefficientGenerator<-function(n,l, mode, randomPeaksBoolean, 
                               randomPeaksFrequency = NULL, randomPeaksRange=NULL, a=NULL,
                               sineWaves = NULL)
{
  #PART I: DEFINE THE REGULAR PART OF THE COEFFICIENT SEQUENCE
  if(mode==1)  #log_2-l decay, i.e. c_j=log2(1+(1/j^(l)))
  {
    C<-matrix(apply(matrix(seq(1:n), nrow=1), 2, function(x){log2(1+(1/x^(l)))}), ncol=1)
  }
  else if(mode==2) #c_j=a^j
  {
    if(missing(a)) #display warning for debugging purposes
    {
      print("WARNING! 'a' was not given as function input in coefficientGenerator even though mode==2!")
    }
    C<-matrix(apply(matrix(seq(1:n), nrow=1), 2, function(x){a^x}))
  }
  else if(mode==0) #c_j=1/j^l decay
  {
    C<-matrix(apply(matrix(seq(1:n), nrow=1), 2, function(x){(1/x^(l))}), ncol=1)
  }
  else if(mode == 3) #sine decay
  {
    #how many waves do we want? sineWaves many.
    multiplier<-sineWaves * (2*pi) #sine is 2-pi periodic
    C<-matrix(apply(matrix(seq(1:n), nrow=1), 2, function(x){sin((x/n)*multiplier)}), ncol=1)
  }
  #PART II: DEFINE THE IRREGULAR PART OF THE COEFFICIENT SEQUENCE
  #If you want an irregular coefficient sequence, random peaks will be injected.
  if(randomPeaksBoolean)
  {
    if(missing(randomPeaksFrequency)|missing(randomPeaksRange))
    {
      print("WARNING! 'randomPeaksFrequency' or 'randomPeaksRange' was not given
            as function input in coefficientGenerator even though randomPeaksBoolean == TRUE!")
    }
    peaks<-runif(n*randomPeaksFrequency, -randomPeaksRange, randomPeaksRange) #generate the peaks
    peakPositions<-sample.int(n,size=n*randomPeaksFrequency)                  #generate the peak positions
    C[peakPositions]<-peaks                                                   #replace respective elements in C
  }
  return(C)
}

#generate a WOLD-type process, return the T observations in X
WoldGenerator<-function(n,T,k,coefSequence)
{
  #Step 1: Generate the errors
  E<-errorGenerator(n,T,k)
  #Step 2: apply the coefficient sequence to the errors and obtain x_t
  X<-matrix(  apply(matrix(seq(1:T), nrow=1), 2, function(x){t(   matrix(E[(x:(x+n-1)),1], ncol=1)   )%*%(coefSequence)}),  ncol=1)
  #Step 3: return the T observations of x_t:0<t<T+1
  return(X)
}

#calculate the true gamma-parameters implied by the generated coefficient sequence, i.e. E(X_tX_{t-h}) using h=1,2,...,p
# Notice that if k implies a t-distribution, the variance is obtained as k/(k-2)
trueGammas<-function(coefSequence, p,k,n)
{
  gamma_0<-t(coefSequence)%*%coefSequence
  gamma_p<-matrix(apply(matrix(seq(1:p), nrow=1), 2, function(x){t(coefSequence[1:(n-x)])%*%coefSequence[(x+1):n]}))
  Gamma_p<-matrix(1,nrow=p,ncol=p)
  sigma2<-1
  for(i in 1:p)
  {
    for(j in 1:p)
    {
      if(i==j)
      {
        Gamma_p[i,j]<-gamma_0
      }
      else
      {
        Gamma_p[i,j]<-gamma_p[abs(i-j)]
      }
    }
  }
  if(k>2)
  {
    sigma2<-(k/(k-2))
  }
  list(gamma_0=gamma_0, gamma_p=gamma_p, Gamma_p = Gamma_p, sigma2=sigma2)
}

trueAlphas<-function(trueGammaList)
{
  alphas<-solve(trueGammaList$Gamma_p)%*%(trueGammaList$gamma_p)
  return(alphas)
}


