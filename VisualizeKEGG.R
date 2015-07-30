setwd("/Users/vivianhsiao/Rprojects/KEGGVis")
source("KEGGio.R")
source("DataUtils.R")
require("Rgraphviz")
require("KEGGgraph")
require("shape")

plotKEGGgraph <- function(pathway.id, organism="mmu", default.color="#99CCFF", node.colors=NULL,
                          node.labels=NULL, save.png=FALSE, save.dir.path="./images/", rank="") {
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
  
  # Get default attributes
  at <- getMyDefaultAttr(default.color="#99CCFF")
  
  # Set node labels
  nAttrs <- list()
  if(!is.null(node.labels)) {
    nAttrs$label <- node.labels
  }
  
  graph <- layoutGraph(graph, nodeAttrs=nAttrs, attrs=at)
  nodeRenderInfo(graph) <- myNodeInfoDefaults(nodes(graph))
  edgeRenderInfo(graph) <- myEdgeInfoDefaults(edges(graph))
  
  # Set custom node colors
  if(!is.null(node.colors)) {
    nodeRenderInfo(graph) <- list(fill=node.colors)
  }
  graphRenderInfo(graph) <- list(main=pathway.name, cex.main=3)
  
  ht = max(length(nodes)*60, 500)
  if(save.png) {
    title <- gsub("\\s+", "_", standardize(pathway.name))
    png(filename = paste0(save.dir.path, "/", rank, organism, pathway.id, title, ".png"),
        height=ht, width=max(1.5*ht, 500), units="px")
  }
  frame()
  renderGraph(graph)
  if(save.png) {
    dev.off()
  }
}

appendToLabels <- function(node.labels, node.label.appends) {
  names <- names(node.labels)
  n.labels <- paste0(node.labels, sep="\n", node.label.appends)
  names(n.labels) <- names
  return(n.labels)
}

myNodeInfoDefaults <- function(nodes) {
  n.i.def <- list(lwd=2)
  return(n.i.def)
}

myEdgeInfoDefaults <- function(edges) {
  e.i.def <- list(arrowhead="vee")
  return(e.i.def)
}

getMyDefaultAttr <- function(default.color = "#99CCFF", default.shape="ellipse") {
  attrs <- list(node=list(shape=default.shape, fixedsize=FALSE, fillcolor=default.color, fontsize=14),
                edge=list(arrowsize=0.2))
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

hot.cold.colors <- function(n, cold="green", hot="red") {
  f <- colorRampPalette(c(cold, "white", hot))
  return(f(n))
}