


#' @title A function to compute Horizontal Visibility Graphs and associated statistics
#' @export
#' @description Calculates a Horizontal Visibility Graph
#' @usage HVG(x, meth, maxL, rho)
#' @param x A time series
#' @param meth A character string that describes the HVG method to use. Currently implemented: "HVG", "HVG_weighted", "LPHVG", "LPHVG_weighted".
#' @param maxL Maximum length of the time series.
#' @param rho Additional parameter
#' @details
#' Horizontal Visibility Graphs map a time series into a complex network.
#' Following Luque, B., Lacasa, L., Ballesteros, F. and Luque, J., 2009. Horizontal visibility graphs: Exact results for random time series. Physical Review E, 80(4), p.046103.
#' ATTENTION: This function is still in development and needs further testing!
#' @return A list that contains the adjacency matrix, degree distribution, and further HVG-based statistics.
#' @references Luque, B., Lacasa, L., Ballesteros, F. and Luque, J., 2009. Physical Review E, 80(4), p.046103.
#' @author Sebastian Sippel
#' @examples
#' x = arima.sim(model=list(ar = 0.3), n = 10^2)
#' HVG(x, meth = "HVG", maxL = 10^9, rho = NA)
HVG <- function(x, meth = "HVG", maxL = 10^9, rho = NA) {
  
  # % Ben Fulcher, October 2009
  # based on C routines
  
  # %% Preliminaries, check inputs
  y <- x; rm(x)  # 
  N = length(y) # ; % time-series length
  
  if ( N > maxL ) {
    message("Time series is too long for visibility graph!")
    return(NULL)
  }
  
  if ( N > 1000 ) {
    # implement sparse matrices here (!)
    A = Matrix::Matrix(data = 0, nrow = N, ncol = N, sparse = T)
    # A_forward = Matrix::Matrix(data = 0, nrow = N, ncol = N, sparse = T)
    # A_backward = Matrix::Matrix(data = 0, nrow = N, ncol = N, sparse = T)
  } else {
    # implement zeros:
    A = matrix(data = 0, nrow = N, ncol = N)
    # A_forward = matrix(data = 0, nrow = N, ncol = N)
    # A_backward = matrix(data = 0, nrow = N, ncol = N)
  }
  
  y = y - min(y)
  
  # Visibility Graphs:
  # ------------------------------------------------
  if (meth == "VG") {
    # for (i in 1:(N-1)) {
    # compute all subsequent gradients:
    #  deltay = 
    
    #}
    message("Not implemented at the moment")
    return(NULL)
  }
  
  # HVG:
  # ------------------------------------------------
  if (meth == "HVG" | meth == "HVG_weighted") {
    rho = 0
    y.block = .C("HVG_C", x = as.double(y), 
                 N = as.integer(N),
                 forwardBlocker = as.integer(rep(-1, N)),
                 backwardBlocker = as.integer(rep(-1, N)))
  }
  
  ## Limited Penetrable HVG:
  # ------------------------------------------------
  if (meth == "LPHVG" | meth == "LPHVG_weighted") {
    y = (y - mean(y)) / stats::sd(y)
    
    y.block = .C("HVG_penetrable_C", x = as.double(y), 
                 N = as.integer(N),
                 forwardBlocker = as.integer(rep(-1, N * (rho+1))),
                 backwardBlocker = as.integer(rep(-1, N * (rho+1))),
                 rho = as.integer(rho))
  }
  
  na.forward = which(y.block$forwardBlocker == -1)
  na.backward = which(y.block$backwardBlocker == -1)
  
  ## For unweighted visibility graph:
  # -----------------------------------------------
  if (meth == "HVG" | meth == "LPHVG") {
    # fill matrix with forward values:
    A[cbind(rep(c(1:N), (rho+1))[-na.forward], y.block$forwardBlocker[-na.forward])] <- 1
    
    # fill matrix with backward values:
    A[cbind(y.block$backwardBlocker[-na.backward], rep(c(1:N), (rho+1))[-na.backward])] <- 1
  }
  
  ## For weighted visibility graph:
  # -----------------------------------------------
  if (meth == "HVG_weighted" | meth == "LPHVG_weighted") {
    # fill matrix with forward values:
    A[cbind(rep(c(1:N), (rho+1))[-na.forward], y.block$forwardBlocker[-na.forward])] <- abs(y[rep(c(1:N), (rho+1))[-na.forward]] - y[y.block$forwardBlocker[-na.forward]])
    
    # fill matrix with backward values:
    A[cbind(y.block$backwardBlocker[-na.backward], rep(c(1:N), (rho+1))[-na.backward])] <- abs(y[rep(c(1:N), (rho+1))[-na.backward]] - y[y.block$backwardBlocker[-na.backward]])
  }
  
  
  # before A is symmetrized:
  k_outgoing = Matrix::rowSums(A)  
  k_incoming = Matrix::rowSums(Matrix::t(A))
  
  # symmetrize A's crudely:
  A_outgoing = A
  A_incoming = Matrix::t(A)
  A <- A + Matrix::t(A)
  
  # Statistics on the output:
  k = Matrix::rowSums(A) # the connectivity time series
  
  if (meth == "HVG" | meth == "LPHVG") {
    kdis = sapply(X = 1:max(k), FUN = function(i) length(which(k == i)) / length(k))
    kdis.outgoing = sapply(X = 1:max(k_outgoing), FUN = function(i) length(which(k_outgoing == i)) / length(k_outgoing))
    kdis.incoming = sapply(X = 1:max(k_incoming), FUN = function(i) length(which(k_incoming == i)) / length(k_incoming))
    
    # get distance distribution as well:
    ddis.forwardBlocker=c(y.block$forwardBlocker - (1:N)) #[-na.forward]
    ddis.forwardBlocker[na.forward] = Inf
    
    ddis.backwardBlocker=c((1:N)-y.block$backwardBlocker)
    ddis.backwardBlocker[na.backward] = Inf
    ddis.all = c(ddis.forwardBlocker, ddis.backwardBlocker)
    dis.idx=sort(unique(ddis.all), decreasing = F)
    Inf.idx=which(ddis.all==Inf)
    
    
    ddis=graphics::hist(x = ddis.all[-Inf.idx], breaks= seq(0.5, to = max(dis.idx[-which(dis.idx==Inf)])+0.5, 1), plot=F)$counts
    ddis.norm=ddis/sum(ddis)
    Inf.ratio=length(Inf.idx)/N
    #ddis1= rep(0, max(ddis.all[-Inf.idx]))
    #ddis1[dis.idx[-which(dis.idx==Inf)]]=sapply(X = dis.idx[-which(dis.idx==Inf)], FUN=function(i) length(which(ddis.all[-Inf.idx] == i)) / length(ddis.all))
    #plot(ddis, log="xy", type='l')
  }
  
  if (meth == "HVG_weighted" | meth == "LPHVG_weighted") {
    kdis = sapply(X = 1:ceiling(max(k)), FUN = function(i) length(which(k < i & k >= (i-1))) / length(k))
    kdis.outgoing = sapply(X = 1:ceiling(max(k_outgoing)), FUN = function(i) length(which(k_outgoing < i & k_outgoing >= (i-1))) / length(k_outgoing))
    kdis.incoming = sapply(X = 1:ceiling(max(k_incoming)), FUN = function(i) length(which(k_incoming < i & k_incoming >= (i-1))) / length(k_incoming))
  }
  
  # get Survival function / complimentary cumulative distribution function:
  ccdf = rev(cumsum(rev(kdis)))
  ccdf.outgoing = rev(cumsum(rev(kdis.outgoing)))
  ccdf.incoming = rev(cumsum(rev(kdis.incoming)))
  
  mode.k = mode(k) 
  mean.k = mean(k)
  median.k = stats::median(k)
  sd.k = stats::sd(k)
  max.k = max(k)
  min.k = min(k)
  range.k =range(k)
  iqr.k = stats::IQR(k)
  max.over.median = max(k) / stats::median(k)
  k.statistics = list(mode.k = mode.k, mean.k = mean.k, 
                      median.k = median.k, sd.k = sd.k,
                      max.k = max.k, min.k = min.k, 
                      range.k = range.k, iqr.k = iqr.k,
                      max.over.median = max.over.median)
  
  return(list(A = A, A.outgoing = A_outgoing, A.incoming = A_incoming, 
              k = k, k.outgoing = k_outgoing, k.incoming = k_incoming,
              kdis = kdis, kdis.outgoing = kdis.outgoing, kdis.incoming = kdis.incoming,
              ccdf = ccdf, ccdf.outgoing = ccdf.outgoing, ccdf.incoming = ccdf.incoming,
              k.statistics = k.statistics, ddis = ddis.norm, Inf.ratio=Inf.ratio ))
}

