
# -----------------------------------------
## Partitioning of time series using symbolic sequences (i.e., median partitioning)...:
# -----------------------------------------

# Sebastian Sippel
# based on: Hauhs and Lange, 2008: http://onlinelibrary.wiley.com/doi/10.1111/j.1749-8198.2007.00075.x/epdf


# Define functions: migfcanalyse
#' @title A function to compute Mean information gain (MIG) and Fluctuation complexity (FC)
#' @export
#' @description Calculates MIG and FC
#' @usage migfc(x, L)
#' @param x A time series
#' @param L word length parameter
#' @details
#' MIG and FC are based on a median partitioning of the time series
#' Following Hauhs, M. and Lange, H., 2008. Classification of runoff in headwater catchments: A physical problem?. Geography Compass, 2(1), pp.235-254.
#' ATTENTION: This function is still in development and needs further testing!
#' @return A list containing MIG, FC and transition matrices.
#' @references Hauhs, M. and Lange, H., (2008) Geography Compass, 2(1), pp.235-254.
#' @author Sebastian Sippel
#' @examples
#' x = arima.sim(model=list(ar = 0.3), n = 10^4)
#' migfc(x, L=4)
migfc <- function(x, L) {  
  
  # get median partitioning:
  symseq = as.numeric(x > stats::median(x))
  # partition into L-words:
  wL=get.sym.words(symseq,L)
  wLp1=get.sym.words(symseq,L+1)
  ret.list = get.MIG_FC(wL=wL,wLp1=wLp1)
  
  return(ret.list)
}

#' @keywords internal
get.sym.words <- function(sym.seq, L) {
  # partition into word counts:
  wolist=rep(0, 2^L)
  
  for (i in 1:(length(sym.seq)-(L-1))) {
    wn=1;
    wo = c()
    for (j in 1:L) {
      wo[j]=sym.seq[i+j-1]
      wn=wn+round(wo[j]*2^(L-j))
    }
    wolist[wn]=wolist[wn]+1
  }
  return(wolist)
}  


# output of "get.sym.words" is input to komp-program:
# wL = get.sym.words(sym.seq= sym.seq, L=3)
# wLp1 = get.sym.words(sym.seq= sym.seq, L=4)

#' @keywords internal
get.MIG_FC <- function(wL, wLp1) {
  n=length(wL)
  if(n != length(wLp1)/2) {
    message('dimensions of wL and wLp1 are not in a 1:2 ratio!')
    return(NULL) }
  
  pijmat = matrix(data=0, nrow = n, ncol=n) # zeros(n,n);pi_jmat=zeros(n,n);
  pi_jmat = matrix(data=0, nrow = n, ncol=n)
  ant=sum(wLp1);
  for (i in (1:n)) {
    for (j in (2*i-1):(2*i)) {
      pijmat[i,j-n*floor(j/(n+0.5))]=wLp1[j]/ant
      pi_jmat[i,j-n*floor(j/(n+0.5))]=wLp1[j]/wL[i]
      # message(j)
    }
  }
  
  # calculate mean information gain and fluctuation complexity
  migi=0;fluc=0;
  for (i in (1:n)) {
    for (j in (1:n)) {
      if(pi_jmat[i,j]>0) {
        ausmig = -pijmat[i,j] * log2(pi_jmat[i,j])   # ; %disp(aus);
        migi = migi+ausmig }
      if(wL[i]>0 && wL[j]>0) {
        ausfc = pijmat[i,j] * log2(wL[i]/wL[j])^2  # ; %disp(aus);
        fluc = fluc+ausfc }
      # message(j)
    }
  }
  ret.list = list(migi, fluc, pijmat, pi_jmat)
  names(ret.list) <- c("MIG", "FC", "pijmat", "pi_jmat")
  return(ret.list)
}  


