standardize.default <- function(call, unchanged = NULL,
                                 standardize.y = FALSE, binary.inputs = "center") {
  form <- call$formula
  varnames <- all.vars(form)
  n.vars <- length(varnames)
  
  # Decide which variables will be unchanged
  transform <- rep("leave.alone", n.vars)
  if (standardize.y) {
    transform[1] <- "full"
  }
  
  for (i in 2:n.vars) {
    v <- varnames[i]
    # Retrieve the data
    if (is.null(call$data)) {
      thedata <- get(v)
    } else {
      thedata <- get(as.character(call$data))[[v]]
    }
    
    # Check if thedata is valid
    if (is.function(thedata)) {
      stop(paste("Error: The object", v, "is a function, not a data vector."))
    }
    if (is.null(thedata) || (!is.numeric(thedata) && !is.factor(thedata))) {
      stop(paste("Error: The object", v, "must be either numeric or a factor."))
    }
    
    if (is.na(match(v, unchanged))) {
      num.categories <- length(unique(thedata[!is.na(thedata)]))
      
      if (num.categories == 2) {
        transform[i] <- binary.inputs  # Handle binary inputs
      } else if (num.categories > 2 & is.numeric(thedata)) {
        transform[i] <- "full"  # Treat numeric with more than 2 categories as standard
      }
    }
  }
  
  # New variable names: prefix with "c." if centered or "z." if scaled
  varnames.new <- ifelse(transform == "leave.alone", varnames,
                          ifelse(transform == "full", paste("z", varnames, sep = "."),
                                 paste("c", varnames, sep = ".")))
  
  transformed.variables <- (1:n.vars)[transform != "leave.alone"]
  
  # Define the new variables
  if (is.null(call$data)) {
    for (i in transformed.variables) {
      assign(varnames.new[i], rescale(get(varnames[i]), binary.inputs))
    }
  } else {
    newvars <- NULL
    for (i in transformed.variables) {
      new_var <- rescale(get(as.character(call$data))[[varnames[i]]], binary.inputs)
      assign(varnames.new[i], new_var)  # Assign to global environment
      newvars <- cbind(newvars, new_var)  # Combine new variables for output
    }
    assign(as.character(call$data), cbind(get(as.character(call$data)), newvars))
  }
  
  # Now call the regression with the new variables
  call.new <- call
  L <- sapply(as.list(varnames.new), as.name)
  names(L) <- varnames
  call.new$formula <- do.call(substitute, list(form, L))
  
  formula <- as.character(call.new$formula)
  if (length(formula) != 3) stop("formula does not have three components")
  
  formula <- paste(formula[2], formula[1], formula[3])
  formula <- gsub("factor(z.", "factor(", formula, fixed = TRUE)
  formula <- gsub("factor(c.", "factor(", formula, fixed = TRUE)
  call.new$formula <- as.formula(formula)
  
  return(eval(call.new))
}


setMethod("standardize", signature(object = "lm"),
  function(object, unchanged = NULL, 
           standardize.y = FALSE, binary.inputs = "center") {
    call <- object$call
    out <- standardize.default(call = call, unchanged = unchanged, 
                               standardize.y = standardize.y, 
                               binary.inputs = binary.inputs)
    return(out)
  }
)

setMethod("standardize", signature(object = "glm"),
  function(object, unchanged = NULL, 
           standardize.y = FALSE, binary.inputs = "center") {
    call <- object$call
    out <- standardize.default(call = call, unchanged = unchanged, 
                               standardize.y = standardize.y, 
                               binary.inputs = binary.inputs)
    return(out)
  }
)

setMethod("standardize", signature(object = "polr"),
  function(object, unchanged = NULL, 
           standardize.y = FALSE, binary.inputs = "center") {
    call <- object$call
    out <- standardize.default(call = call, unchanged = unchanged, 
                               standardize.y = standardize.y, 
                               binary.inputs = binary.inputs)
    return(out)
  }
)

setMethod("standardize", signature(object = "merMod"),
  function(object, unchanged = NULL, 
           standardize.y = FALSE, binary.inputs = "center") {
    call <- object@call  # For merMod, use @ instead of $
    out <- standardize.default(call = call, unchanged = unchanged, 
                               standardize.y = standardize.y, 
                               binary.inputs = binary.inputs)
    return(out)
  }
)
