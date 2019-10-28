## Sebastian Sippel
# 11.03.2018

# This script is to write functions that calculate 
# complexity measures based on the Renyi entropy
# see M. Jauregui et al., Physica A, 498 74-85, 2018

# ----------------------------------------------------
# Permutation Entropy based on Renyi:
# ----------------------------------------------------
# test = rnorm(10^4)
# test_opd = ordinal_pattern_distribution(x = test, ndemb = 4)


#' @title A function to compute Renyi entropy
#' @export
#' @description Renyi permutation entropy
#' @usage permutation_entropy_Renyi(opd, alpha)
#' @param opd A numeric vector that details an ordinal pattern distribution.
#' @param alpha alpha parameter in Renyi entropy
#' @details
#' This function calculates the Renyi entropy as described in Jauregui et al., Physica A, 498 74-85, 2018.
#' @return The Renyi entropy value.
#' @references Jauregui et al., Physica A, 498 74-85, 2018.
#' @author Sebastian Sippel
#' @examples
#' x = arima.sim(model=list(ar = 0.3), n = 10^4)
#' opd = ordinal_pattern_distribution(x = x, ndemb = 6)
#' permutation_entropy_Renyi(opd = opd, alpha = 0.5)
permutation_entropy_Renyi = function(opd, alpha) {
  if (alpha <= 0) return(NULL)
  if (alpha == 1) return(permutation_entropy(opd = opd))
  
  # maximum Shannon Entropy is uniform distribution:
  ssmax  = log( length(opd) )
  # compute Shannon entropy based on Renyi:
  S_alpha = shannon_entropy_Renyi(opd = opd, alpha = alpha)
  H_alpha = S_alpha / ssmax
  return(H_alpha)
}



# compute non-normalized Shannon Entropy:
#' @keywords internal
shannon_entropy_Renyi = function(opd, alpha) {
  opd.prob = opd / sum(opd)
  S_alpha = 1 / (1 - alpha) * log(sum(sapply(opd.prob, FUN=function(prob) if (prob >= 1.e-30) return(prob ^ alpha) else return(0))))
  return(S_alpha)
}



# ----------------------------------------------------
# q_complexity:
# ----------------------------------------------------

#' @title A function to compute Renyi complexity
#' @export
#' @description Renyi complexity
#' @usage complexity_Renyi(opd, alpha)
#' @param opd A numeric vector that details an ordinal pattern distribution.
#' @param alpha alpha parameter in Renyi complexity
#' @details
#' This function calculates the Renyi complexity as described in Jauregui et al., Physica A, 498 74-85, 2018.
#' @return The Renyi complexity value.
#' @references Jauregui et al., Physica A, 498 74-85, 2018.
#' @author Sebastian Sippel
#' @examples
#' x = arima.sim(model=list(ar = 0.3), n = 10^4)
#' opd = ordinal_pattern_distribution(x = x, ndemb = 6)
#' complexity_Renyi(opd = opd, alpha = 0.5)
complexity_Renyi = function(opd, alpha) {
  
  
  # convert to probabilities:
  opd.prob = opd/sum(opd)
  opd.length = length(opd)
  
  if (alpha == 1) return(MPR_complexity(opd = opd))
  
  # H_q (Normalized Shannon entropy  based on q-log):
  H_alpha = permutation_entropy_Renyi(opd = opd, alpha = alpha)
  
  
  D_alpha = 1 / (2 * (alpha - 1)) * (
    log(sum(sapply(opd.prob, FUN=function(prob) prob ^ alpha * ((prob + 1 / opd.length) / 2) ^ (1-alpha) ))) + 
      log(sum(sapply(opd.prob, FUN=function(prob) 1 / opd.length ^ alpha * ((prob + 1 /opd.length) / 2) ^ (1-alpha)))) )
  
  D_alpha_max = 1 / ( 2 * (alpha-1)) * log( ((opd.length + 1)^(1-alpha) + opd.length - 1 ) / opd.length * ((opd.length + 1)/(4*opd.length)) ^ (1-alpha))
  
  C_alpha = D_alpha * H_alpha / D_alpha_max
  
  return(C_alpha)
}

