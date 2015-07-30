setwd("/Users/vivianhsiao/Rprojects/KEGGVis")

combineAnalyses <- function(analysis.files="example", analysis.types="example", by.column=1,
                            join='outer') {
  # Return an outer join of analyses in analysis.files.
  tryCatch({
    if(analysis.files == "example") {
      analysis.files <- c("./example_data/example_analysis.csv", "./example_data/example_fc.csv")
    }
    if(analysis.types == "example") {
      analysis.types <- c("t.test", "fold.change")
    }
  }, warning = function(w) {})

  if(length(analysis.files < 1)) {
    simpleError("No analysis files were given.")
  } 
  tryCatch({
    first = TRUE
    for (file in analysis.files) {
      newf <- read.table(file, header=TRUE, sep=",", quote="", colClasses="character")
      if(first) {
        combined <- newf
        first = FALSE
      } else {
        combined <- merge(combined, newf, all=TRUE, by = by.column, suffixes="_2")
      }
    }
  })
  return(combined)
}