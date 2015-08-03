setwd("/Users/vivianhsiao/Rprojects/KEGGVis")

combineAnalyses <- function(analysis.files, analysis.types, by.column=1,
                            join='outer') {
  # Return an outer join of analyses in analysis.files.

  if(length(analysis.files < 1)) {
    simpleError("No analysis files were given.")
  } 
  tryCatch({
    first = TRUE
    for (file in analysis.files) {
      newf <- read.table(file, header=TRUE, sep=",", colClasses="character")
      if(first) {
        combined <- newf
        first = FALSE
      } else {
        combined <- merge(combined, newf, all=TRUE, by = by.column)
      }
    }
  })
  return(combined)
}