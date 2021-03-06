#' compute the optimal changepoint model for a vector of real-valued data
#' and a non-negative real-valued penalty,
#' given the square loss (to minimize) / gaussian likelihood (to maximize)
#'
#' @param data A numerical vector for which the changepoint model is to be computed
#' @param penalty A non-negative real number indicating penalty parameter
#' @return A vector of the optimal cost values and a vector of the optimal segment ends
#' @examples
#' data(neuroblastoma, package="neuroblastoma")
#' selectedData <- subset(neuroblastoma$profiles, profile.id=="1" & chromosome=="1")
#' opart::opart_gaussian(selectedData$logratio, 1)
#' @export


opart_gaussian <- function(data, penalty) {

  #check for positive length of data
  if(!(length(data) > 0)){
    stop("data vector must have atleast one element")
  }

  #check if there are any missing values in data vector
  if(any(is.na(data))){
    stop("data vector has missing(NA) values")
  }

  #check if data vector has all finite numeric values
  if(!(all(is.numeric(data)) && all(is.finite(data)))){
    stop("data vector must contains finite numeric values")
  }

  #check if the penalty value is numeric
  if(!is.numeric(penalty)){
    stop("penalty value should be numeric")
  }

  #check if penalty is finite and of length 1
  if(!(is.finite(penalty) && (length(penalty) == 1))){
    stop("penalty must be a finite numeric value")
  }

  result <- .C("opart_gaussian_interface",
               n_data = as.integer(length(data)),
               data.vec = as.double(data),
               penalty = as.double(penalty),
               cost.vec = double(length(data)),
               sums = double(length(data)),
               dp = double(length(data)),
               end.vec = integer(length(data)),
               positions = integer(length(data)),
               PACKAGE="opart")

  seg_ends <- (result$end.vec)

  #remove -2 placeholders from the output
  result$end.vec <- seg_ends[seg_ends != -2]

  #remove the columns used for internal calculations as they don't need to be displayed
  result <- result[c("cost.vec","end.vec")]

  #display the result
  result
}
