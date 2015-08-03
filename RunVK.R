setwd("/Users/vivianhsiao/Rprojects/KEGGVis")
source("VisualizeKEGG.R")
require("Rgraphviz")
require("KEGGgraph")
require("shape")

renderAffectedPathwaysCombined <- function(pathway.list.file="example", organism="mmu", alias.file="./cached.pathways/globalpathways.alias.csv", 
                                           analysis.paths="example", analysis.types=c("t.test", "Fold.Change"), color.according.to="Fold.Change", 
                                           spectrum=hot.cold.colors(100), label.cols=c("p.value", "Fold.Change")) {
  # Set up example files
  example.dir <- "./example_data"
  if(analysis.paths == "example") {
    analysis.paths = sapply(c("/example_analysis.csv", "/example_fc.csv"), function (x) paste0(example.dir, x))
  }
  if(pathway.list.file == "example") {
    pathway.list.file <- paste0(example.dir, "/example_pathway_list.csv")
  }
  
  analysis.dir <- dirname(analysis.paths[1])
  combined <- combineAnalyses(analysis.files=analysis.paths, analysis.types=analysis.types)
  combined.path <- paste0(analysis.dir, "/combined.csv")
  write.table(combined, file=combined.path, sep = ",", row.names = FALSE)
  renderAffectedPathways(pathway.list.file = pathway.list.file, organism = organism, alias.file = alias.file, 
                         analysis.path = combined.path, analysis.type = "combined", color.according.to = color.according.to,
                         spectrum = spectrum, label.cols = label.cols)
}

renderAffectedPathways <- function(pathway.list.file="example", organism="mmu", alias.file="./cached.pathways/globalpathways.alias.csv", 
                                   analysis.path="example", analysis.type="t.test", color.according.to="p.value",
                                   spectrum=heat.colors(100), label.cols=c("p.value")) {
  
  print("Rendering affected pathways...")
  # Set up example files
  example.dir <- "./example_data"
  if(analysis.path == "example") {
    analysis.path = paste0(example.dir, "/example_analysis.csv")
  }
  if(pathway.list.file == "example") {
    pathway.list.file <- paste0(example.dir, "/example_pathway_list.csv")
  }
  
  alias.file.path = paste0("./cached.pathways/", alias.file)
  pathways <- vector()
  tryCatch({
    pathways <- readLines(pathway.list.file)
  }, error = function(e) {
    simpleMessage(paste("Problem reading pathway list file -", pathway.list.file))
    simpleMessage(e)
  })
  
  # Fetch aliases from file
  al <- read.table(alias.file, header=TRUE, sep=",", strip.white = TRUE, colClasses = "character")
  al$aliases2 <- strsplit(al$aliases, ";")
  al$aliases2 <- sapply(al$aliases2, function(x) standardize(x))
  
  pathway.ids <- list()
  rank <- 1
  ranks <- list()
  
  for(pway in pathways) {
    pway <- standardize(pway)
    alias.id = find.alias(pway, al, rank)
    if(!is.null(alias.id)) {
      ranks <- c(ranks, rank)
      pathway.ids <- c(pathway.ids, alias.id)
    } else {
      print(sprintf("No alias found for %s", pway))
    }
    rank <- rank + 1;
  }
  for(j in 1:length(pathway.ids)) {
    id <- pathway.ids[j]
    rank <- ranks[j]
    
    # Set up directories
    save.dir.path <- paste0(dirname(analysis.path), "/images")
    if(!dir.exists(save.dir.path)) {
      dir.create(save.dir.path)
    }
    
    renderAffectedPathway(id, organism=organism, analysis.path = analysis.path, analysis.type = analysis.type,
                          color.according.to = color.according.to, spectrum = spectrum, save.png = TRUE,
                          save.dir.path=save.dir.path, rank=rank, label.cols = label.cols)
  }
}

renderAffectedPathway <- function(pathway.id="00020", organism="mmu", analysis.path="./example_data/example_analysis.csv", 
                                  analysis.type="t.test", color.according.to="p.value", spectrum=heat.colors(100),
                                  save.png=FALSE, save.dir.path="./images/", rank="", label.cols = c("p.value")) {
  print(paste("Rendering:", pathway.id))
  analysis <- data.frame()
  compounds <- vector()
  if(analysis.type %in% c("t.test", "combined")){
    analysis <- read.table(analysis.path, sep=",", row.names = 1, header=TRUE)
    compounds <- paste0("cpd:", row.names(analysis))
    p.value <- analysis$p.value
  } else {
    simpleError(paste("Analysis type not supported:", analysis.type))
  }
  
  # Color according to data
  node.colors <- color.data(analysis, color.according.to, spectrum)
  node.labels <- label.data(pathway.id, organism=organism, analysis, label.cols)
  
  plotKEGGgraph(pathway.id, organism=organism, node.colors = node.colors, node.labels = node.labels,
                save.png = save.png, save.dir.path = save.dir.path, rank=rank)
}

# Helper functions
find.alias <- function(pathway, al, rank) {
  alias.id <- NULL
  for(i in 1:nrow(al)) {
    id <- al$id[[i]]
    aliases <- al$aliases2[[i]]
    if(pathway %in% aliases) {
      print(sprintf("Alias found for %s: %s [rank:%d]", pathway, aliases[1], rank))
      alias.id <- id
      break()
    } 
  }
  return(alias.id)
}

label.data <- function(pathway.id, organism=organism, analysis, label.cols, with.id=FALSE) {
  node.labels = list()
  names <- row.names(analysis)
  # Set labels to be values in name map
  tryCatch({
    name.map = getAndCacheNameMap(pathway.id, organism=organism)
    node.labels = as.character(name.map$V2)
    node.labels = formatLabels(node.labels)
    names = as.character(name.map$V1)
    names(node.labels) <- names
  }, error = function(e) {
    print(paste("ERROR:", e))
  })
  
  # Add data to labels
  l.data <- get.label.data(analysis, label.cols, with.id=with.id)
  for (l.d.n in names(l.data)) {
    cpd <- paste0("cpd:",l.d.n)
    if(cpd %in% names(node.labels)) {
      node.labels[cpd] <- paste0(node.labels[cpd], "\n\n", l.data[l.d.n])
    }
  }
  return(node.labels)  
}

get.label.data <- function(analysis, label.cols, with.id=FALSE, alpha = 0.05, show.Holm = TRUE) {
  tryCatch({
    first = TRUE
    for (dc in label.cols) {
      dat = unlist(analysis[dc])
      if (dc=="p.value") {
        # Holm P
        holm.p = p.adjust(dat, method = "holm")
        dc.data = paste0(dc, ": ", sprintf("%0.2e (Holm: %0.2e)", dat, holm.p))
        
        # Find significant entries
        sig.idx <- (dat < alpha & !is.na(dat))
        sig.idx.hp <- (holm.p < alpha & !is.na(holm.p))
        
        # Add * for significant, ** if Holm adjusted p is significant
        dc.data[sig.idx] = paste("*", dc.data[sig.idx])
        dc.data[sig.idx.hp] = paste0("*", dc.data[sig.idx.hp])
      } else {
        dc.data = paste0(dc, ": ", sprintf("%0.2f", dat))
      }
      if(first) {
        l.data <- dc.data
        first = FALSE
      } else {
        l.data <- paste0(l.data, sep="\n", dc.data)
      }
    }
    if(with.id) {
      l.data <- paste0(l.data, sep="\n", row.names(analysis))
    }
    names(l.data) <- row.names(analysis)
  }, warning = function(w) {
    simpleWarning(paste("Could not add parameter to label", label.cols, "- default labels only."))
    l.data = list()
  }, error = function(e) {
    simpleError(e)
  })
  return(l.data)
}

color.data <- function(analysis, color.according.to, spectrum){
  tryCatch({
    color.data <- unlist(analysis[color.according.to])
    if(color.according.to %in% c("p.value")) {
      color.data <- -log10(color.data)
    } else if(color.according.to %in% c("Fold.Change")) {
      color.data <- log2(color.data)
    }
    
    if (color.according.to %in% c("Fold.Change")) {
      max.abs <- max(abs(color.data), na.rm = TRUE)
      mid <- ceiling(length(spectrum) / 2)
      idx <- ceiling(mid + ((color.data / max.abs) * mid))
      node.colors = spectrum[idx]
    } else {
      max.color <- max(color.data)
      node.colors = spectrum[ceiling((color.data / max.color) * length(spectrum))]
    }
  }, warning = function(w) {
    simpleWarning(paste("Could not color according to parameter", color.according.to, "- default coloring only."))
    node.colors = NULL
  })
  names(node.colors) <- paste0("cpd:", row.names(analysis))
  return(node.colors)
}


formatLabels <- function(labels, n=15) {
  sapply(labels, formatLabel)
}

formatLabel <- function(label, n=20) {
  label <- strsplit(label, ";")[[1]]
  label <- label[nchar(label) == min(nchar(label))]
  label <- breakLineAtNCharsAndSpaces(label)
  return(label)
}

breakLineAtNCharsAndSpaces <- function(line, n=30) {
  words <- strsplit(strip(line), "\\s", perl=TRUE)[[1]]
  buf <- vector()
  for (w in words) {
    for (i in seq(1, nchar(w), n)) {
      buf <- c(buf, substr(w, i, i+n-1))
    }
  }
  new <- paste0(buf, collapse="\n")
  return(new)
}

strip <- function(str) {
  stripped <- gsub("^\\s+|\\s+$", "", str, perl = TRUE)
  return(stripped)
}

standardize <- function(x) {
  std <- x
  std <- tolower(std)
  std <- gsub("[[:punct:]]", "", std) # Remove punctuation
  std <- gsub("\\s+", " ", std) # remove spaces >1 long
  std <- strip(std) # Remove leading/trailing ws
  return(std)
}