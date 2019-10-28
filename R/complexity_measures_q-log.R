## Sebastian Sippel
# 11.03.2018

# This script is to write functions that calculate 
# complexity measures based on the q-logarithm
# see H. V. Ribeiro et al., 2017, ArXiv "Characterizing Time Series via Complexity-Entropy Curves
# "
# https://arxiv.org/abs/1705.04779

# ----------------------------------------------------
# Permutation Entropy based on q-log:
# ----------------------------------------------------

#' @title A function to compute q-log permutation entropy
#' @export
#' @description q-log permutation entropy
#' @usage permutation_entropy_qlog(opd, q)
#' @param opd A numeric vector that details an ordinal pattern distribution.
#' @param q q-log parameter
#' @details
#' This function calculates the q-log permutation entropy as described in Ribeiro et al. 2017.
#' @return The q-log permutation entropy value.
#' @references Ribeiro et al. 2017, https://arxiv.org/abs/1705.04779.
#' @author Sebastian Sippel
#' @examples
#' x = arima.sim(model=list(ar = 0.3), n = 10^4)
#' opd = ordinal_pattern_distribution(x = x, ndemb = 6)
#' permutation_entropy_qlog(opd = opd, q = 1)
permutation_entropy_qlog = function(opd, q) {
  # maximum Shannon Entropy is uniform distribution:
  # ssmax  = log( length(opd) )
  ssmax = q_log(x = length(opd), q = q)

  # compute Permutation entropy and return:
  PE_qlog = shannon_entropy_qlog(opd = opd, q = q) / ssmax 
  return(PE_qlog)
}

# compute non-normalized Shannon Entropy:
#' @keywords internal
shannon_entropy_qlog = function(opd, q) {
  opd.prob = opd / sum(opd)
  H_s = sum(sapply(opd.prob, FUN=function(prob) if (prob >= 1.e-30) return(prob * (q_log(x = 1 / prob, q = q))) else return(0)))
  return(H_s)
}
# test = rnorm(10000)
# test_opd = ordinal_pattern_distribution(test, ndemb = 5)
# shannon_entropy_qlog(opd = test_opd, q = 1)


# get q-logarithm
#' @keywords internal
q_log <- function(x, q) {
  if ( x == 0) return(0)
  
  if (q != 1) {
    return((x^(1-q)-1) / (1-q))
  } else if (q == 1) {
    return(log(x))
  }
}



# ----------------------------------------------------
# q_complexity:
# ----------------------------------------------------
#' @title A function to compute q-log complexity
#' @export
#' @description q-log complexity
#' @usage q_complexity(opd, q)
#' @param opd A numeric vector that details an ordinal pattern distribution.
#' @param q q-log parameter
#' @details
#' This function calculates the q-log complexity as described in Ribeiro et al. 2017.
#' @return The q-log complexity value.
#' @references Ribeiro et al. 2017, https://arxiv.org/abs/1705.04779.
#' @author Sebastian Sippel
#' @examples
#' x = arima.sim(model=list(ar = 0.3), n = 10^4)
#' opd = ordinal_pattern_distribution(x = x, ndemb = 6)
#' q_complexity(opd = opd, q = 1)
q_complexity = function(opd, q) {
  
 
  # convert to probabilities:
  opd.prob = opd/sum(opd)
  opd.length = length(opd)
  
  if (q == 1) return(MPR_complexity(opd = opd))
  
  # H_q (Normalized Shannon entropy  based on q-log):
  H_q = permutation_entropy_qlog(opd = opd, q = q)
  
  # D_q: (according to Ribeiro et al. Phys. A; Divergence measure between P and U):
  opd.prob.check = opd.prob[which(opd.prob >= 1.e-30)]
  
  D_q_P_U = - 0.5 * sum(sapply(X = opd.prob.check, FUN = function(p_i) p_i * q_log(x = (p_i + 1 / opd.length) / ( 2 * p_i ), q = q))) -
    0.5 * sum(sapply(X = opd.prob, FUN = function(p_i) 1 / opd.length * q_log(x = (p_i + 1 / opd.length) / ( 2 / opd.length), q = q ) ))
  
  D_q_max = ( (2 ^ (2 - q)) * opd.length - (1+opd.length) ^ (1-q) - opd.length * (1 + 1 / opd.length)^(1-q) - opd.length + 1 ) / ((1-q) * 2^(2-q) * opd.length)
  
  C_q = D_q_P_U * H_q / D_q_max
  
  return(C_q)
}
