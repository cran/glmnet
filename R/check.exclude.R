check.exclude <- function(exclude,nvars){
      exclude <- unique(exclude)
      if(length(exclude)> (nvars-2))stop("cannot retain 1 or less variables")
      if(length(exclude)==0)exclude <- NULL
      exclude
      }
