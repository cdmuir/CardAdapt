# R functions associated with Muir and Angert 2016

# Make statistical significance asterisks
sigStar <- function(pvalue)
{
  if (pvalue >= 0.05) return("")
  if (pvalue < 0.05 & pvalue >= 0.01) return("*")
  if (pvalue < 0.01 & pvalue >= 0.001) return("**")
  if (pvalue < 0.001) return("***")
}

# Determine a number's order of magnitude
oom <- function(x)
{
  # x should be a vector of numbers
  floor(log10(abs(x)))
}

# Determine number of significant digits to use for rounding in tables
sigDig <- function(x)
{
  # x should be a vector of numeric elements to be rounded to same significant digit
  ret <- round(x, -min(oom(x)))
  return(ret)
}

# Text string for p-value in tables
pval2latex <- function(x)
{
  stem <- x * 10 ^ -oom(x)
  ret <- ifelse(oom(x) > -3, 
                sprintf("%.*f", 3, x), 
                sprintf("%.*f $\\\\times10^{%s}$", 1, stem, oom(x)))
  return(ret)
}


