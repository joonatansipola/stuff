## Tools

# Get indexes of matching rows from 2 matrices
match.rows <- function(x, y, all=F) 
{
  if (all) 
  {
    if(ncol(x) != ncol(y)) return(rep(F, nrow(x)))
    ret <- list()
    for(x.row in 1:nrow(x)) {
      x.row.matches <- vector()
      for(y.row in 1:nrow(y)) {
        if(all(x[x.row,] == y[y.row,])) {
          x.row.matches <- append(x.row.matches, y.row)
        }
      }
      ret[[x.row]] <- x.row.matches
    }
    return(ret)
  }
  
  else {
    x <- apply(x, 1, paste, collapse="cOlLaPsE")
    y <- apply(y, 1, paste, collapse="cOlLaPsE")
    match(x,y)
  }
}

# Grep vector with multiple patterns and return indexes
double.grep.idx <- function(pattern1, pattern2, search.from) 
{
  which(grepl(search.from, pattern=pattern1) & grepl(search.from, pattern=pattern2))
}
triple.grep.idx <- function(pattern1, pattern2, pattern3, search.from) 
{
  which(grepl(search.from, pattern=pattern1) & grepl(search.from, pattern=pattern2) & grepl(search.from, pattern=pattern3))
}

# return vector without NAs or NaNs
na.remove <- function(vec) {
  vec[!is.na(vec) & !is.nan(vec) & vec != "NaN"]
}