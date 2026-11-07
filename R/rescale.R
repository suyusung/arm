rescale <- function(x, binary.inputs = "center") {
    # Convert x to numeric if its not already
    if (!is.numeric(x)) {
        # Store the original levels for categorical variables
        levels_x <- levels(factor(x))
        x <- as.numeric(factor(x))
        
        # Handle binary factors by mapping levels to 0 and 1
        if (length(levels_x) == 2) {
            x <- x - 1  # Maps two levels to 0 and 1
        }
    }
    
    x.obs <- x[!is.na(x)]
    
    if (length(unique(x.obs)) == 2) {
        # For binary factors
        if (binary.inputs == "0/1") {
            x <- (x-min(x.obs))/(max(x.obs)-min(x.obs))
            return(x)  # Return original scale
        }
        else if (binary.inputs == "-0.5,0.5") {
            # Map the two observed values onto -0.5 and 0.5, whatever
            # they happen to be, so a variable coded 1/2 is treated
            # like one coded 0/1.  Subtracting 0.5 from x only works
            # when the input is already 0/1, which gave 0.5/1.5 for a
            # 1/2-coded variable (the same bug as issue #19, which was
            # fixed for factors only).
            return((x - min(x.obs)) / (max(x.obs) - min(x.obs)) - 0.5)
        }
        else if (binary.inputs == "center") {
            return(x - mean(x.obs))  # Center the variable
        }
        else if (binary.inputs == "full") {
            return((x - mean(x.obs)) / (2 * sd(x.obs)))  # Standardize
        }
    } else {
             return ((x-mean(x.obs))/(2*sd(x.obs)))
        }
}

# Example usage
