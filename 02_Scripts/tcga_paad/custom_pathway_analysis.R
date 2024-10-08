scriptLibraries <-  c(
  "DESeq2",
  "ggplot2",
  "ggrepel",
  "RColorBrewer",
  "DT",
  "pheatmap",
  "limma", 
  "tidyverse",
  "here",
  "ggpubr",
  "sjPlot",
  "Biobase",
  "clusterProfiler",
  "org.Hs.eg.db",
  "vroom",
  "msigdbr",
  "DOSE",
  "ReactomePA",
  "GOstats",
  "genefilter",
  "topGO",
  "fgsea",
  "GOplot",
  "gprofiler2",
  "igraph",
  "intergraph",
  "GGally",
  "tm",
  "wordcloud",
  "piano",
  "biomaRt"
)

# scripts --------------------------------------------------------------
devtools::source_url("https://raw.githubusercontent.com/umahajanatlmu/useful_commands/main/basicFunctions.R")
# load packages -------------------------------------------------------
installScriptLibs(scriptLibraries)



##----------------------------------------------------------------
##                       topGo functions                       --
##----------------------------------------------------------------
myTopGoBP <- function(foreG, res, species = c("hsapiens", "mmusculus")) {
  
  if (species == "hsapiens") {
    map <- "org.Hs.eg.db"
  } else if (species == "mmusculus") {
    map <- "org.Mm.eg.db"
  }
  
  foreG <- na.omit(foreG)
  res <- res %>%
    drop_na(EnsemblID)
  ## row names to entrez
  overallBaseMean <- as.matrix(res[, "baseMean", drop = F])
  sig_idx <- match(foreG, res$EnsemblID)
  backG <- c()
  for(i in sig_idx){
    ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
    backG <- c(backG, ind)
    
  }
  backG <- unique(backG)
  backG <- res$EnsemblID[backG]
  
  backG <- setdiff(backG,  foreG)
  
  geneIDs <- res$EnsemblID
  inUniverse <- geneIDs %in% c(foreG,  backG) 
  inSelection <-  geneIDs %in% foreG
  alg <- factor( as.integer(inSelection[inUniverse]))
  #alg <- factor( as.integer( inSelection ) )
  names(alg) <- geneIDs[inUniverse]
  #names(alg) <- geneIDs
  
  tgd <- new( "topGOdata", ontology="BP", allGenes = alg, nodeSize=5,
              annot=annFUN.org, mapping=map, ID = "ensembl" )
  
  ## run tests
  resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
  resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
  
  ## look at results
  bp <- GenTable( tgd, Fisher.elim = resultTopGO.elim, 
                  Fisher.classic = resultTopGO.classic,
                  orderBy = "Fisher.classic" , topNodes = 100)
  return(bp)
  
}

myTopGoMF <- function(foreG, res, species = c("hsapiens", "mmusculus")) {
  
  if (species == "hsapiens") {
    map <- "org.Hs.eg.db"
  } else if (species == "mmusculus") {
    map <- "org.Mm.eg.db"
  }
  
  foreG <- na.omit(foreG)
  res <- res %>%
    drop_na(EnsemblID)
  overallBaseMean <- as.matrix(res[, "baseMean", drop = F])
  sig_idx <- match(foreG, res$EnsemblID)
  backG <- c()
  for(i in sig_idx){
    ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
    backG <- c(backG, ind)
    
  }
  backG <- unique(backG)
  backG <- res$EnsemblID[backG]
  
  backG <- setdiff(backG,  foreG)
  
  geneIDs = res$EnsemblID
  inUniverse = geneIDs %in% c(foreG,  backG) 
  inSelection =  geneIDs %in% foreG
  alg <- factor( as.integer(inSelection[inUniverse]))
  #alg <- factor( as.integer( inSelection ) )
  names(alg) <- geneIDs[inUniverse]
  #names(alg) <- geneIDs
  
  
  tgd <- new( "topGOdata", ontology="MF", allGenes = alg, nodeSize=5,
              annot=annFUN.org, mapping=map, ID = "ensembl" )
  
  ## run tests
  resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
  resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
  
  ## look at results
  mf <- GenTable( tgd, Fisher.elim = resultTopGO.elim, 
                  Fisher.classic = resultTopGO.classic,
                  orderBy = "Fisher.classic" , topNodes = 100)
  return(mf)
}

myTopGoCC <- function(foreG, res, species = c("hsapiens", "mmusculus")) {
  
  if (species == "hsapiens") {
    map <- "org.Hs.eg.db"
  } else if (species == "mmusculus") {
    map <- "org.Mm.eg.db"
  }
  
  foreG <- na.omit(foreG)
  res <- res %>%
    drop_na(EnsemblID)
  overallBaseMean <- as.matrix(res[, "baseMean", drop = F])
  sig_idx <- match(foreG, res$EnsemblID)
  backG <- c()
  for(i in sig_idx){
    ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
    backG <- c(backG, ind)
    
  }
  backG <- unique(backG)
  backG <- res$EnsemblID[backG]
  
  backG <- setdiff(backG,  foreG)
  
  geneIDs = res$EnsemblID
  inUniverse = geneIDs %in% c(foreG,  backG) 
  inSelection =  geneIDs %in% foreG
  alg <- factor(as.integer(inSelection[inUniverse]))
  #alg <- factor( as.integer( inSelection ) )
  names(alg) <- geneIDs[inUniverse]
  #names(alg) <- geneIDs
  
  tgd <- new( "topGOdata", ontology="CC", allGenes = alg, nodeSize=5,
              annot=annFUN.org, mapping=map, ID = "ensembl" )
  
  ## run tests
  resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
  resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
  
  ## look at results
  cc <- GenTable( tgd, Fisher.elim = resultTopGO.elim, 
                  Fisher.classic = resultTopGO.classic,
                  orderBy = "Fisher.classic" , topNodes = 100)
  return(cc)
}

##----------------------------------------------------------------
##                     keggStats functions                     --
##----------------------------------------------------------------
keggStats <- function(foreground, background, hgCutoff = 0.001, species =c("hsapiens", "mmusculus")) {
  
  if (species == "hsapiens") {
    annotation = "org.Hs.eg.db"
    header = "hsa"
  } else if (species == "mmusculus") {
    annotation = "org.Mm.eg.db"
    header = "mmu"
  }
  
  entrezUniverse <- na.omit(background)
  selectedEntrezIds <- na.omit(foreground)
  params <- new(
    "KEGGHyperGParams",
    geneIds = selectedEntrezIds,
    universeGeneIds = entrezUniverse,
    annotation = annotation,
    pvalueCutoff = hgCutoff,
    testDirection = "over"
  )
  table <- hyperGTest(params)
  table <- summary(table)
  
  grouped_rows <- split(seq_len(nrow(table)), rep(1:ceiling(nrow(table)/10), each = 10, length.out = nrow(table)))
  
  for (i in grouped_rows) {
    table[i, "Pathway"] <- sapply(KEGGREST::keggGet(paste0(header,table[i, "KEGGID"])), "[[", "NAME")
  }
  
  return(table)
}

##----------------------------------------------------------------
##                     GSEAStats functions                     --
##----------------------------------------------------------------
GSEAStats <- function(df, species =c("hsapiens", "mmusculus"), 
                      GSEA_terms = c("H","c2", "c3", "c4", "c5", "c6", "c7"), 
                      top_n = 5) {
  
  banner(GSEA_terms)
  if (species == "hsapiens") {
    load(url(paste0("https://bioinf.wehi.edu.au/software/MSigDB/human_",GSEA_terms,"_v5p2.rdata")))
    assign("gsea_set",get(paste0("Hs.",GSEA_terms)), envir = .GlobalEnv) 
  } else if (species == "mmusculus") {
    load(url(paste0("https://bioinf.wehi.edu.au/software/MSigDB/mouse_",GSEA_terms,"_v5p2.rdata")))
    assign("gsea_set",get(paste0("Mm.",GSEA_terms)), envir = .GlobalEnv) 
  }
  
  ranks <- df$log2FoldChange
  names(ranks) <- df$EntrezID
  ranks <- ranks[!duplicated(names(ranks))]
  
  fgseaRes <- fgsea(gsea_set, ranks, minSize = 5, maxSize = 500, nperm = 1000)
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n = top_n), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n = top_n), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  
  p <- plotGseaTable(gsea_set[topPathways], ranks, fgseaRes, gseaParam = 0.5)
  
  ## plot gsea
  p_dat <- fgseaRes %>%
    as.data.frame() %>%
    arrange(desc(NES)) %>%
    filter(padj < 0.05) %>%
    filter(row_number() > max(row_number()) - 5 | row_number() <= 5) %>%
    mutate(NES.col = ifelse(NES > 0, "#E41A1C", "#377EB8")) %>%
    mutate(pathway = gsub("_", " ", pathway))
  
  p1 <- ggplot(p_dat, aes(x = reorder(pathway, NES), y = NES)) + 
    geom_col(aes(fill = NES.col),
             color = "black", 
             show.legend = FALSE) + 
    coord_flip() + 
    labs(x = "",
         y = "Normalized Enrichment Score", 
         title = paste("GSEA:",GSEA_terms)) +
    theme_classic()+
    theme(plot.title = element_text(face = "bold", size = 12),
          legend.position = "bottom", 
          axis.text.x = element_text(size = 12), 
          axis.text.y = element_blank(), 
          axis.title = element_text(face = "bold",
                                    size = 12),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank()
    ) + 
    scale_fill_identity() +
    geom_text(
      aes(label = reorder(pathway, NES),
          y = ifelse(p_dat$NES < 0, 0.1, -0.1)),
      position = position_dodge(width = 0.9),
      hjust = ifelse(p_dat$NES < 0, 0, 1),
      color = "black",
      size = 3
    ) +
    geom_hline(yintercept = 0)
  
  result <- list(gsea_p = p,
                 plot = p1,
                 result=topPathways)
  
  return(result)
  
}

##----------------------------------------------------------------
##                 network annotation function                 --
##----------------------------------------------------------------
## ref: https://github.com/tfkillian/GSEA_network_analysis/blob/master/GSEA_network_analysis.Rmd
t.rW <-
  c(
    "cell",
    "process",
    "regulation",
    "negative",
    "positive",
    "signaling",
    "response",
    "stimulus",
    "signal",
    "activity",
    "protein",
    "involved",
    "component",
    "level",
    "effector",
    "event",
    "projection",
    "organismal",
    "cellular",
    "modification",
    "pathway",
    "mediated",
    "dependent",
    "organization",
    "group",
    "target",
    "biocarta",
    "kegg",
    "reactome"
  )
clust_head <- function(x){
  txt <- unlist(strsplit(x, "_"))
  txt <- Corpus(VectorSource(txt))
  txt <- tm_map(txt, PlainTextDocument)
  txt <- tm_map(txt, removePunctuation)
  txt <- tm_map(txt, removeNumbers)
  txt <- tm_map(txt, content_transformer(tolower))
  txt <- tm_map(txt, removeWords, c(t.rW, stopwords("english")))
  tdm <- TermDocumentMatrix(txt)
  m <- as.matrix(tdm)
  word_freqs <- sort(rowSums(m), decreasing=TRUE)
  word_freqs <- word_freqs[word_freqs>1]
  word_freqs <- paste(names(word_freqs)[1:4], collapse=" ")
  gsub("[[:space:]]?NA[[:space:]]?", "", word_freqs)
}

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

##----------------------------------------------------------------
##                       network function                       --
##----------------------------------------------------------------

networkStats <- function(df, species =c("hsapiens", "mmusculus"), 
                         top_n = 3,
                         p_val = 0.05) {
  terms <- c()
  for (i in c("H","c2", "c3", "c4", "c5", "c6", "c7")) {
    if (species == "hsapiens") {
      load(url(paste0("https://bioinf.wehi.edu.au/software/MSigDB/human_",i,"_v5p2.rdata")))
      assign("gsea_set",get(paste0("Hs.",i)), envir = .GlobalEnv) 
    } else if (species == "mmusculus") {
      load(url(paste0("https://bioinf.wehi.edu.au/software/MSigDB/mouse_",i,"_v5p2.rdata")))
      assign("gsea_set",get(paste0("Mm.",i)), envir = .GlobalEnv) 
    }
    terms <- append(terms, gsea_set)
  }
  gs.libs <- sapply(names(terms), function(x) strsplit(x, "_")[[1]][1])
  gset <- terms[which(gs.libs %in% c("KEGG", "REACTOME", "BIOCARTA", "GO", "HALLMARK"))]
  
  idx <- ids2indices(gene.sets = gset, identifiers = df$EntrezID)
  
  dat <- cameraPR(df$log2FoldChange, idx, sort = FALSE)
  dat$PValue.Mixed <- cameraPR(abs(df$log2FoldChange), idx, sort = FALSE)$PValue
  dat$FDR.Mixed <- p.adjust(dat$PValue.Mixed, method = "BH")
  dat$name <- rownames(dat)
  
  dat$Direction <- as.character(dat$Direction)
  dat$Direction[dat$FDR > p_val] <- "Mixed"
  dat$Direction[dat$Direction == "Mixed" & dat$FDR.Mixed > p_val] = "NOT"
  dat$Direction <- factor(dat$Direction, levels = c("NOT", "Up", "Down", "Mixed"))
  
  idx <- which(dat$Direction == "Mixed")
  if (length(idx) > 0)
    dat$FDR[idx] = dat$FDR.Mixed[idx]
  
  dat <- dat[, -grep("\\.Mixed", names(dat))]
  dat <- dat[dat$Direction != "NOT", ]
  dat$Direction <- factor(dat$Direction, levels = c("Up", "Down", "Mixed"))
  
  # only keep gene sets present in the data
  id.keep <- which(names(gset) %in% dat$name)
  gset.subset <- gset[id.keep]
  
  # adjacency matrix
  m.adj <- sapply(gset.subset, 
                  function(x) 
                    sapply(gset.subset, function(y) length(intersect(unlist(x),unlist(y)))))
  
  diag(m.adj) = 0
  
  # Jaccard index matrix
  NGenes <- sapply(gset.subset, length)
  m.union <- outer(NGenes, NGenes, "+") - m.adj
  m.jacc <- m.adj/m.union
  
  # apply cutoff to Jaccard matrix
  m.adj1 <- m.adj * (m.jacc > p_val)
  
  # construct network object
  net <- graph_from_adjacency_matrix(m.adj1, "max", diag = FALSE, weighted = TRUE)
  
  # add vertex features
  V(net)$size <- dat$NGenes
  
  # choose node colors
  palette <- c("#e41a1c","#1b9e77","#569FC9")
  names(palette) <- c("Up", "Down", "Mixed")
  
  V(net)$color <- palette[dat$Direction]
  V(net)$Direction <- as.character(dat$Direction)
  
  singletons <- which(igraph::degree(net) == 0)
  net1 <- delete_vertices(net, singletons)
  in.single <- which(dat$name %in% V(net)$name[singletons])
  tab.net <- dat[in.single, ]
  tab.net$FDR <- signif(tab.net$FDR, 2)
  tab.net$name <- gsub("_", " ", tab.net$name)
  
  ## identify binary systems
  
  clu1 <- igraph::components(net1)
  clu.lt3 <- which(sizes(clu1) < top_n)
  v.clu.lt3 <- which(clu1$membership %in% clu.lt3)
  net2 <- delete_vertices(net1, v.clu.lt3)
  clu2 <- igraph::components(net2)
  in.clu.lt3 <- which(dat$name %in% V(net1)$name[v.clu.lt3])
  tab.net1 <- dat[in.clu.lt3, ]
  tab.net1$FDR <- signif(tab.net1$FDR, 2)
  cludp <- clu1$membership[v.clu.lt3]
  cludp <- data.frame(name = names(cludp), id = as.numeric(cludp))
  tab.net1 <- merge(tab.net1, cludp)
  tab.net1$name <- gsub("_", " ", tab.net1$name)
  
  ## detect communities
  net2 <- delete_edge_attr(net2, "weight")
  clu3 <- cluster_edge_betweenness(net2)
  
  # delete edges between communities
  net3 <- delete_edges(net2, which(as.vector(crossing(clu3, net2))))
  
  # remove clusters of size <5
  small_cluster_ids <- which(sizes(clu3) < top_n)
  small_cl_v <- which(clu3$membership %in% small_cluster_ids)
  net3 <- delete_vertices(net3, small_cl_v)
  
  clu3 <- igraph::components(net3)
  
  nodecol <- c(brewer.pal(9, "Set1"), brewer.pal(9, "Dark2"))
  nodecol <- colorRampPalette(nodecol)(max(clu3$membership))
  
  ## lattice of annotated networks
  clust <- data.frame(cl = clu3$membership)
  rownames(clust) <- names(V(net3))
  
  # generate cluster titles
  cl3.lab.txt <- tapply(rownames(clust), clust$cl, clust_head)
  
  # remove NAs
  cl3.lab.txt <- gsub("[[:space:]]?NA[[:space:]]?", "", cl3.lab.txt)
  cl3.lab.txt <- data.frame(cl3.lab.txt)
  # clu3 <- igraph::components(net3)
  clu.order <- order(clu3$csize, decreasing = TRUE)
  clu3$mem <- match(clu3$membership, clu.order)
  
  clu3$label <- cl3.lab.txt$cl3.lab.txt[clu3$membership]
  clu3$label[duplicated(clu3$label)] <- NA
  clu3$label <- toupper(clu3$label)
  clu3$label <- sub("(^.{10,30})[[:space:]]", "\\1\\\n", clu3$label)
  
  ## print
  p <- ggnet2(net3, 
              size = 0, 
              label = FALSE, 
              color = nodecol[clu3$membership],
              edge.size = 1, 
              edge.color = "grey") + 
    geom_point(color = "black") +
    geom_point(aes(color = color)) + 
    geom_label_repel(aes(label = clu3$label,
                         color = color), 
                     alpha = 0.75, 
                     fontface = "bold", 
                     size = 3, 
                     box.padding = 0.1,
                     label.padding = 0.1, 
                     na.rm = TRUE)
  
  # generate a list of ggplots
  g <- list(max(clu3$membership))
  set.seed(123456)
  
  for (ii in 1:max(clu3$membership)) {
    subgf <- induced_subgraph(net3, which(clu3$mem == ii))
    
    # generate titles with one optional line break
    title <- substr(toupper(cl3.lab.txt$cl3.lab.txt[clu.order][ii]), 1, 60)
    
    if (nchar(title) > 25) {
      title <- sub("(^.{10,30})[[:space:]]", "\\1\\\n", title)
    }
    
    # generate node labels using word 2-5 of the geneset name
    
    v.label <- names(V(subgf))
    
    v.label <- lapply(v.label, function(x) strsplit(x, "_")[[1]])
    v.label <- sapply(v.label, function(x) paste(x[2:min(5, length(x))], collapse = "_"))
    
    # clean up geneset names
    v.label <- gsub("_PATHWAY", "", v.label)
    v.label <- gsub("_SIGNALING", "", v.label)
    
    # introduce line breaks
    v.label <- gsub("_", "\n", v.label)
    
    # remove node labels for large clusters
    if (length(v.label) > 5)
      v.label = rep(NA, length(v.label))
    
    g[[ii]] = ggnet2(subgf, 
                     edge.size = 0.5, 
                     edge.color = "#66c2a5", label = FALSE,
                     size = V(subgf)$size, 
                     max_size = 3, 
                     size.cut = 4, 
                     color = palette[V(subgf)$Direction]) +
      theme(legend.position = "none", 
            plot.title = element_text(size = 6),
            panel.grid = element_blank()) + 
      geom_label_repel(label = v.label,
                       size = 2, 
                       box.padding = 0.1, 
                       label.padding = 0.1) + 
      ggtitle(title)
    
  }
  
  nr.cols <- min(4, max(clu3$membership))
  nr.rows <- ceiling(max(clu3$membership)/nr.cols)
  
  width <- sapply(g, function(x) nrow(x$data))
  grid.arrange <- getFromNamespace("grid.arrange", asNamespace("gridExtra"))
  
  if (length(g) < 20) {
    g.length <- length(g)
  } else 
    g.length <- 20
  
  g <- grid.arrange(grobs = g[seq(g.length)], 
                    ncol = nr.cols)
  
  return(list(map = p,
              ind = g))
  
}
