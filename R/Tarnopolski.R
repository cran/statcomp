

## Tarnopolski diagrams and turning point:
# ------------------------------------------------
# Tarnopolski diagram: plot Abbe-value vs. Turning point ratio
# following M. Tarnopolski / Physica A 461 (2016) 662-673
# Sebastian Sippel, 18.10.2019

# install.packages("somebm")
# require(somebm)

# 2.2 Ratio of mean square successive difference to the variance
# following M. Tarnopolski / Physica A 461 (2016) 662-673
# so-called "Abbe"-value, quantifies "the degree of smoothness of a time series"
# decreases to zero for very smooth time series, tends to unity for purely white noise

#' @title A function to compute Abbe values
#' @export
#' @description Calculates "Abbe" values.
#' @usage Abbe(x)
#' @param x A time series
#' @details
#' "Abbe" values quantify the degree of smoothness of a time series. 
#' Decreases to zero for very smooth time series, tends to unity for purely white noise.
#' Following Tarnopolski et al. Physica A 461 (2016) 662-673.
#' @return Abbe value
#' @references Tarnopolski et al. (2016), Physica A 461, 662-673.
#' @author Sebastian Sippel
#' @examples
#' x = arima.sim(model=list(ar = 0.3), n = 10^4)
#' Abbe(x)
Abbe <- function(x) {
  n=length(x)
  d2=1/(n-1) * sum(diff(x)^2)
  s2=1/n * sum((x-mean(x))^2)
  Abbe = d2 / (2 * s2)
  return(Abbe)
}


# 2.3 Turning point analysis:
# x=rnorm(10)
# plot(x, type='b')

#' @title A function to compute Turning points
#' @export
#' @description Calculates Turning point values.
#' @usage Turning_point(x)
#' @param x A time series
#' @details
#' Turning point of a time series.
#' Following Tarnopolski et al. Physica A 461 (2016) 662-673.
#' @return Turning point value
#' @references Tarnopolski et al. (2016), Physica A 461, 662-673.
#' @author Sebastian Sippel
#' @examples
#' x = arima.sim(model=list(ar = 0.3), n = 10^4)
#' Turning_point(x)
Turning_point <- function(x) {
  n=length(x)
  TuP = length(which(abs(diff(sign(diff(x))))==2))
  E_TuP = 2/3 * (n - 2)
  return(TuP / E_TuP)
}


# Turning_point(x = cumsum(rnorm(100000)))
# plot Abbe-value vs. Turning point ratio for Tarnopolski diagram:

# fbm=lapply(X = seq(0.1, 0.5, 0.1), FUN=function(x) fbm(hurst = x, n = 10^4))
# plot(c(1,1), type='n', xlab = "Abbe", ylab ="Turning point rate", xlim = c(0, 1.5), ylim = c(0.7, 1.2))
# lapply(X = fbm, FUN=function(x) points(Abbe(x), Turning_point(x), xlim = c(0, 1.5), cex = 1, pch = 16, col = "darkblue"))
# x=rnorm(10^4)
# points(x=Abbe(x), y = Turning_point(x), pch = 16)
# x=runif(10^4)
# points(x=Abbe(x), y = Turning_point(x), pch = 16)

