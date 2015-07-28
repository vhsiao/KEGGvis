setwd("/Users/vivianhsiao/Rprojects/KEGGVis")
require("Rgraphviz")
require("KEGGgraph")

getGlobalMap <- function(organism = "mmu", pathway.list.file = 'globalpathways', clear.cache = FALSE, clear.kgml.caches = FALSE) {
  global.map.path <- paste0("./cached.pathways/dot/", pathway.list.file, '.dot')
  map.list.path <- paste0("./cached.pathways/", pathway.list.file) # For now, just grab from file (TODO: grab from online)
  if (clear.cache | !file.exists(global.map.path)) {
    entries <- readLines(file(map.list.path))
    graphs <- list()
    for (ent in entries) {
      if (startsWithKEGGid(ent)) {
        print(paste("Current pathway category:", ent))
      } else {
        parts <- strSplitOnNSpaces(ent, 2)
        pathway.id <- parts[1]
        pathway.label <- parts[2]
        print(paste("id:", pathway.id, "label:", pathway.label))
        
        # Fetch the pathway from KEGG
        tryCatch({
          kgml <- getAndCacheKGML(pathway.id = pathway.id, organism = organism, clear.cache = clear.kgml.caches)
          pathway <- parseKGML(kgml)
          #n <- nodes(kgml)
          tryCatch({
            g <- KEGGpathway2reactionGraph(kgml)  # returns a KEGGgraph
            graphs = c(graphs, list(g))
            }, warning = function(w) {
              print(w)
            })
         }, error = function(e) {
           print(paste("ERROR:", e, "...skipping"))
           traceback(e)
           stop()
         })
      }
    }
    global.graph <- mergeGraphs(graphs, edgemode = "directed")
    # Save to file
    toDot(global.graph, global.map.path)
  } else {
    global.graph <- agread(global.map.path, layoutType = 'dot', layout = TRUE)
  }
  return(global.graph)
}

makeAll <- function(getCacheFunction, organism="mmu", pathway.list.file="globalpathways", clear.cache = FALSE) {
  # Create dot files with edge and node labels
  # Returns a list of pathway id's whose names are pathway names.
  map.list.path <- paste0("./cached.pathways/", pathway.list.file)
  entries <- readLines(file(map.list.path))
  
  # For alias file
  id = vector()
  aliases = vector()
  for (ent in entries) {
    if (startsWithKEGGid(ent)) {
      print(paste("Current pathway category:", ent))
    } else {
      parts <- strSplitOnNSpaces(ent, 2)
      pathway.id <- parts[1]
      pathway.label <- parts[2]
      
      # Add to alias file
      id = c(id, pathway.id)
      aliases = c(aliases, pathway.label)
      print(paste("id:", pathway.id, "label:", pathway.label))
      tryCatch({
        if (!is.null(getCacheFunction)){
          getCacheFunction(pathway.id = pathway.id, organism = organism, clear.cache = FALSE)
        }
        }, error = function(e) {
          print(paste("ERROR:", e))
        })
    }
  }
  if(is.null(getCacheFunction)) {
    # Write alias file
    alias.file = paste0(map.list.path, ".alias.csv")
    alias <- data.frame(id, aliases)
    write.table(alias, file=alias.file, quote=TRUE, sep=",", row.names = FALSE)
  }
  closeAllConnections()
}

resetAliasFile <- function(organism="mmu", pathway.list.file="globalpathways", clear.cache=FALSE) {
  # Generate alias file as tsv
  makeAll(NULL, organism=organism, pathway.list.file=pathway.list.file, clear.cache = clear.cache)
}

makeDotFiles <- function(organism="mmu", pathway.list.file="globalpathways", clear.cache = FALSE) {
  makeAll(getAndCacheDot, organism=organism, pathway.list.file=pathway.list.file, clear.cache = clear.cache)
}

makeNameMaps <- function(organism="mmu", pathway.list.file="globalpathways", clear.cache = FALSE) {
  makeAll(getAndCacheNameMap, organism=organism, pathway.list.file=pathway.list.file, clear.cache = clear.cache)
}

getAndCacheNameMap <- function(pathway.id, organism="mmu", clear.cache = FALSE) {
  destfile = getCachePath(pathway.id, organism="mmu", type="name.map")
  if(clear.cache || !file.exists(destfile)) {
    # Request name info from KEGG
    tryCatch({
      graph <- getKEGGGraph(pathway.id, organism="mmu")
      n <- nodes(graph)
      }, warning = function(w) {
        print(paste0("WARNING:", w))
      })
    # Eg: http://rest.kegg.jp/list/cpd:C00674+cpd:C00117
    url = "http://rest.kegg.jp/list/"
    url = paste0(url, paste0(n, collapse="+"))
    
    # Download the file at the url (will be tsv)
    download.file(url, destfile, "wget")
  }
  # Return data frame of name map
  df <- read.table(destfile, header=FALSE, sep="\t", quote = "")
  return(df)
}

getAndCacheDot <- function(pathway.id, organism="mmu", clear.cache = FALSE) {
  # Convert a KEGG pathway to a dot file with labeled edges and nodes.
  #destfile = paste0("./cached.pathways/dot/", organism, pathway.id, '.dot')
  destfile = getCachePath(pathway.id, organism="mmu", type="dot")
  if(clear.cache || !file.exists(destfile)) {
    pathway <- getAndCacheKGML(pathway.id=pathway.id, organism=organism, clear.cache = FALSE)
    pathway <- parseKGML(pathway)
    g <- KEGGpathway2reactionGraph(pathway)
    if(!is.null(g)) {
      toDot(g, destfile)
      print(paste("Cached dot file for pathway:", pathway.id))
    } else {
      print(paste("No chemical reactions in", pathway.id))
    }
  }
  tryCatch({
    g <- agread(destfile, layoutType = 'dot', layout = TRUE)
  }, warning = function(w) {
    print(paste("Warning:", w))
  }, error = function(e) {
    print(paste("ERROR:", e))
    g <- NULL
  })
  return(g)
}

getKEGGGraph <- function(pathway.id, organism="mmu") {
  graph <- getAndCacheKGML(pathway.id=pathway.id, organism=organism)
  graph <- parseKGML(graph)
  graph <- KEGGpathway2reactionGraph(graph)
  return(graph)
}

getAndCacheKGML <- function(pathway.id, organism="mmu", clear.cache = FALSE) {
  #destfile = paste0("./cached.pathways/xml/", organism, pathway.id, '.xml')
  destfile = getCachePath(pathway.id, organism=organism, type="kgml")
  if(clear.cache || !file.exists(destfile)) {
    pathway <- retrieveKGML(pathwayid = paste0(organism,pathway.id), organism = organism, destfile=destfile)
  }
  return(destfile)
}

getTestGraph <- function(pathway.id = "00030", organism="mmu", pathway=FALSE, reaction=TRUE) {
  tmp <- tempfile()
  g <- retrieveKGML(pathwayid = paste0(organism,pathway.id), organism = organism, destfile=tmp)
  g <- parseKGML(g)
  if(pathway) return(g)
  if(reaction) {
    g <- KEGGpathway2reactionGraph(g)  # returns a KEGGgraph
  } else {
    g <- KEGGpathway2Graph(g, genesOnly=FALSE, expandGenes=FALSE)
  }
  return(g)
}

# Helper Functions

getCachePath <- function(pathway.id, organism="mmu", type="kgml") {
  base.path = "./cached.pathways/"
  paths = list(
    "kgml"=paste0(base.path, "xml/", organism, pathway.id, ".xml"),
    "dot"=paste0(base.path, "dot/", organism, pathway.id, '.dot'),
    "name.map"=paste0(base.path, "name_map/", organism, pathway.id, ".tsv")
  )
  return(paths[type][[1]])
}

startsWithKEGGid <- function(ent) {
  return(regexpr("^[0-9]{5}", ent, perl = TRUE) < 1)
}

strSplitOnNSpaces <- function(ent, n) {
  return(strsplit(ent, sprintf("[\\s]{%d}", n), perl = TRUE)[[1]])
}