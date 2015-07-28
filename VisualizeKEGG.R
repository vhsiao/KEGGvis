setwd("/Users/vivianhsiao/Rprojects/KEGGVis")
source("KEGGio.R")
require("Rgraphviz")
require("KEGGgraph")
require("shape")

renderAffectedPathways <- function(pathway.list.file="example", organism="mmu", alias.file="./cached.pathways/globalpathways.alias.csv", 
                                   analysis.path="example", analysis.type="t.test", color.according.to="p.value",
                                   spectrum=heat.colors(100)) {
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
  al$aliases2 <- sapply(strsplit(al$aliases, ";"), function(x) standardize(x))
  
  pathway.ids <- list()
  rank <- 1
  ranks <- list()
  
  for(pway in pathways) {
    pway <- standardize(pway)
    alias.found = FALSE
    for(i in 1:nrow(al)) {
      id <- al$id[[i]]
      aliases <- al$aliases2[[i]]
      if(pway %in% aliases) {
        print(sprintf("Alias found for %s: %s [rank:%d]", pway, aliases[1], rank))
        pathway.ids <- c(pathway.ids, id)
        alias.found = TRUE
        break()
      } 
    }
    if(alias.found) {
      ranks <- c(ranks, rank)
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
                          save.dir.path = save.dir.path, rank=rank)
  }
}

renderAffectedPathway <- function(pathway.id="00020", organism="mmu", analysis.path="./example_data/example_analysis.csv", 
                                  analysis.type="t.test", color.according.to="p.value", spectrum=heat.colors(100),
                                  save.png=FALSE, save.dir.path="./images/", rank="") {
  analysis <- data.frame()
  compounds <- vector()
  if(analysis.type == "t.test"){
    analysis <- read.table(analysis.path, sep=",", row.names = 1, header=TRUE)
    compounds <- paste0("cpd:", row.names(analysis))
    p.value <- analysis$p.value
  } else {
    simpleError(paste("Analysis type not supported:", analysis.type))
  }
  
  # Color according to data
  tryCatch({
    if(color.according.to == "p.value") {
      color.data <- -log10(p.value)
    } else {
      color.data <- analysis[color.according.to]
    }
    max.color <- max(color.data)
    node.colors = spectrum[ceiling((color.data / max.color) * length(spectrum))]
    #node.colors = rep(spectrum[length(spectrum)], length(compounds))
    names(node.colors) = compounds
  }, warning = function(w) {
    simpleWarning(paste("Could not color according to parameter", color.according.to, "- default coloring only."))
    node.colors = NULL
  })
  
  # Append analysis data to node labels
  # TODO
  node.label.appends = NULL
  plotKEGGgraph(pathway.id, organism=organism, node.colors = node.colors, node.label.appends = node.label.appends,
                save.png = save.png, save.dir.path = save.dir.path, rank=rank)
}

plotKEGGgraph <- function(pathway.id, organism="mmu", default.color="#99CCFF", node.colors=NULL,
                          node.label.appends=NULL, save.png=FALSE, save.dir.path="./images/", rank="") {
  # Get the reaction graph
  tryCatch({
    graph <- getAndCacheKGML(pathway.id=pathway.id, organism=organism)
    pathway <- parseKGML(graph)
    graph <- KEGGpathway2reactionGraph(pathway)
    }, error = function(e) {
      print(paste("ERROR:", e))
  })
  
  pathway.name = getTitle(pathway)
  nodes = nodes(graph)
  
  # Set labels to be values in name map
  tryCatch({
    name.map = getAndCacheNameMap(pathway.id, organism=organism)
    node.labels = as.character(name.map$V2) #TODO: break text if exceeds line width
    node.labels = formatLabels(node.labels)
    names(node.labels) = as.character(name.map$V1)
  }, error = function(e) {
    print(paste("ERROR:", e))
    node.labels = list()
  })
  
  # Set node labels
  nAttrs <- list(label=node.labels)
  # Get default attributes
  at <- getMyDefaultAttr(default.color="#99CCFF")
  
  # Append information to labels
  if(!is.null(node.label.appends)) {
    nAttrs$label <- appendToLabels(node.labels, node.label.appens)
  }
  
  graph <- layoutGraph(graph, nodeAttrs=nAttrs, attrs=at)
  nodeRenderInfo(graph) <- myNodeInfoDefaults(nodes(graph))
  
  # Set custom node colors
  if(!is.null(node.colors)) {
    nodeRenderInfo(graph) <- list(fill=node.colors)
  }
  graphRenderInfo(graph) <- list(main=pathway.name)
  
  ht = max(length(nodes)*50, 500)
  if(save.png) {
    title <- gsub("\\s+", "_", standardize(pathway.name))
    png(filename = paste0(save.dir.path, "/", rank, organism, pathway.id, title, ".png"),
        height=ht, width=max(0.75*ht, 400), units="px")
  }
  frame()
  renderGraph(graph)
  if(save.png) {
    dev.off()
  }
}

myNodeInfoDefaults <- function(nodes) {
  n.i.def <- list(lwd=0.01)
  return(n.i.def)
}
getMyDefaultAttr <- function(default.color = "#99CCFF", default.shape="ellipse") {
  attrs <- list(node=list(shape=default.shape, fixedsize=FALSE, fillcolor=default.color, fontsize=16))
  return(attrs)
}

dispMyDefaultAttr <- function(V=letters[1:10], M=1:4, p=0.2) {
  g <- randomGraph(V, M, p)
  attrs <- getMyDefaultAttrs()
  frame()
  plot(g, attrs=attrs)
}

getTestGraph <- function() {
  return(randomGraph(V=letters[1:10], M=1:4, p=0.2))
}

# Helper functions
appendToLabels <- function(original.labels, appended.info) {
  #TODO
  return(original.labels)
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

standardize <- function(x) {
  std <- x
  std <- tolower(std)
  std <- gsub("[[:punct:]]", "", std) # Remove punctuation
  std <- gsub("\\s+", " ", std) # remove spaces >1 long
  std <- strip(std) # Remove leading/trailing ws
  return(std)
}

strip <- function(str) {
  stripped <- gsub("^\\s+|\\s+$", "", str, perl = TRUE)
  return(stripped)
}