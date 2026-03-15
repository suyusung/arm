rescale <- function(x, binary.inputs = "center") {
    # Convert x to numeric if it's not already
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
            return(x)  # Return original scale
        }
        else if (binary.inputs == "-0.5,0.5") {
            return(x - 0.5)  # Rescale to [-0.5, 0.5]
        }
        else if (binary.inputs == "center") {
            return(x - mean(x.obs))  # Center the variable
        }
        else if (binary.inputs == "full") {
            return((x - mean(x.obs)) / (2 * sd(x.obs)))  # Standardize
        }
    } else {
        # For multinomial variables
        if (binary.inputs == "center") {
            return(x - mean(x.obs))  # Center around mean
        } else if (binary.inputs == "full") {
            return((x - mean(x.obs)) / (2 * sd(x.obs)))  # Standardize
        } else {
            # Default logic for non-binary cases (rescale to [0, 1])
            min_val <- min(x.obs)
            max_val <- max(x.obs)
            return((x - min_val) / (max_val - min_val))  # Rescale to [0, 1]
        }
    }
}

# Example usage
