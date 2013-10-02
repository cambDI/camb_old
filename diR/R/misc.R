# writes a dataframe to a file with a header 
write.with.header <- function(x, file, header, f = write.csv, ...) {
  datafile <- file(file, open = 'wt')
  on.exit(close(datafile))
  if(!missing(header)) writeLines(header,con=datafile)
  f(x, datafile,...)
}


# does dim and head at the same time but with a head of 3
dh <- function(d) {
  print(dim(d))
  print(d[1:3,])
}