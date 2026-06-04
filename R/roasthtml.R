#' Execute ROAST Gene Set Analysis and Output Nested HTML Reports
#'
#' Evaluates multiple collections of gene sets across multiple contrasts using
#' parallelized ROAST operations, generating structured downstream HTML summaries.
#'
#' @param x A list of lists containing `roastgsa` results objects structured by contrast and gene set collection.
#' @param y A list of data.frames containing Differential Expression (DE) metrics structured by contrast.
#' @param mat A numeric matrix of normalized gene expression values.
#' @param mc.cores.x Integer. Number of cores for processing contrasts. Default 1.
#' @param mc.cores.y Integer. Number of cores for processing gene set collections. Default 1.
#' @param outdir Character. Output path for results directories. Default './'.
#' @param maxgs Integer. Maximum number of top paths to include per report. Default 50.
#' @param indhtml Logical. If TRUE, creates dedicated nested child HTML sheets for individual genes. Default TRUE.
#' @param DEdir Character. Source directory path containing static figures (e.g., stripcharts). Default NULL.
#' @param returnData Logical. If TRUE, returns the nested results list structure. Default TRUE.
#' @param intvar Character. Name of the primary design variable of interest. Default NULL.
#' @param mycol Vector of color definitions for heatmap generation.
#'
#' @export
roastHtmlTables <- function(x,y,mat,mc.cores.x=1,mc.cores.y=1,outdir='./',maxgs=50,indhtml=TRUE,DEdir=NULL,returnData=TRUE, intvar=NULL, mycol=mycol)
{
  # check names x == names y
  # check DEdir exists
  ans <- mclapply(names(x[[1]]),function(gs) ## Genesets
  {
    ans <- mclapply(names(y), function(i) ## Contrasts
    {
      mygs <- x[[i]][[gs]]
      ans <- roastHtmlTable(mygs,gs,i,sprintf('roastGSA_MaxMean_%s_%s.html',gs,i),file.path(outdir,sprintf('roastGSA/html/%s/',gs)), indhtml=indhtml,DEdir=DEdir,detable=y[[i]],maxgs=maxgs, mat=mat, intvar=intvar, selcols=colnames(y[[i]]), vorder=colnames(y[[i]])[ncol(y[[i]])], mycol=mycol)
      ans
    },mc.cores=mc.cores.x)
  },mc.cores=mc.cores.y)
  if (returnData) return(ans)
}

#' Construct a Single Gene Set Collection Level HTML Overview
#'
#' Filters gene set statistics using significance thresholds, orchestrates diagnostic
#' image plotting, and calls high-level HTML assembly functions.
#'
#' @param mygs A discrete `roastgsa` collection output list containing `$res` and `$index`.
#' @param gs Character. Target name/identifier of the gene set collection.
#' @param i Character. Target name/identifier of the contrast.
#' @param filename Character. Target HTML filename output.
#' @param dirname Character. Core target path where files will be compiled.
#' @param indhtml Logical. If TRUE, creates individual nested child HTML sheets for genes.
#' @param DEdir Character. Source directory path containing static figures. Default NULL.
#' @param detable A data.frame containing DE table metrics for contrast `i`. Default NULL.
#' @param maxgs Integer. Maximum number of paths to include after filtering. Default 50.
#' @param verbose Logical. If TRUE, logs generation progress to the console. Default TRUE.
#' @param mat A numeric matrix of normalized gene expression values.
#' @param intvar Character. Design variable of interest. Default NULL.
#' @param selcols Character vector. Explicit columns to keep from `detable`.
#' @param vorder Character. Column name to use for descending absolute sorting.
#' @param apval.cut Numeric. Significance threshold value for FDR adjustment filtering. Default 0.05.
#' @param mycol Vector of color definitions for heatmap generation.
#'
#' @export
roastHtmlTable <- function(mygs,gs,i,filename,dirname,indhtml,DEdir=NULL,detable=NULL,maxgs=50,verbose=TRUE, mat, intvar=NULL, selcols, vorder, apval.cut=0.05, mycol=mycol)
{
  if (indhtml & is.null(detable)) stop('indhtml is TRUE but not DE table provided')
  ## Filter for significance, if more than maxgs, maxgs
  mygs$res <- mygs$res[order(abs(mygs$res$nes),decreasing=TRUE),]
  sel <- mygs$res$adj.pval < apval.cut
  if(sum(sel) < 10){
    mmax <- min(10,nrow(mygs$res))
    sel <- 1:mmax
  }
  mygs$res <- mygs$res[sel, ]
  mygs$res <- mygs$res[1:min(nrow(mygs$res),maxgs),]
  mygs$index <- mygs$index[rownames(mygs$res)]
  if (verbose) print(sprintf('Writing output for %s in %s (%d selected genesets)',gs,i,nrow(mygs$res)))
  ## Write down main pathway table file
  if (indhtml) geneDEhtmlfiles <- sprintf('%s_indhtml/%s_genes.html',i,rownames(mygs$res)) else geneDEhtmlfiles <- NULL
  ##dir.create(file.path(dirname,sprintf('roastGSA/html/%s',gs)),recursive=TRUE)
  dir.create(dirname,recursive=TRUE)
  htmlrgsa2(mygs,htmlname=sprintf('roastGSA_MaxMean_%s_%s.html',gs,i),htmlpath=file.path(tablesdir,sprintf('roastGSA/html/%s/',gs)),
            plotpath=sprintf('%s_images/',i),indheatmap=FALSE,y=mat,intvar=intvar,ploteffsize=FALSE,mycol=mycol,
            geneDEhtmlfiles=geneDEhtmlfiles,sorttable=sorttable,dragtable=dragtable,whplot=rownames(mygs$res)[1:nrow(mygs$res)],
            title=sprintf('<center><h4>%s | %s</h4></center>',gs,i))
  if (indhtml) roastHtmlDETable(mygs,gs,i,detable,DEdir,filename,dirname,maxgs,verbose, selcols=selcols, vorder=vorder)
}

#' Create Individual Gene-Level DE Detail HTML Sheets
#'
#' Extracts matching genes for a specific pathway, rounds numeric metrics, sets up
#' JavaScript sorting/dragging wrappers, and compiles individual child HTML sheets.
#'
#' @param mygs A discrete `roastgsa` collection output list.
#' @param gs Character. Name of the gene set collection.
#' @param i Character. Name of the contrast.
#' @param detable A data.frame containing DE table metrics for contrast `i`.
#' @param DEdir Character. Path to base directory containing stripcharts.
#' @param filename Character. Main index filename string.
#' @param dirname Character. Core directory path where nested files will be created.
#' @param maxgs Integer. Maximum number of paths to process. Default 50.
#' @param verbose Logical. If TRUE, logs internal loop progress. Default TRUE.
#' @param selcols Character vector. Explicit columns to subset from `detable`.
#' @param vorder Character. Numeric column name to sort gene rows by.
#'
#' @export
roastHtmlDETable <- function(mygs,gs,i,detable,DEdir,filename,dirname,maxgs=50,verbose=TRUE, selcols, vorder)
{
  outdir <- file.path(dirname,sprintf('%s_indhtml',i))
  dir.create(outdir,recursive=TRUE)
  st <- system.file("javascript", "sorttable.js", package = "phenoTest")
  dt <- system.file("javascript", "dragtable.js", package = "phenoTest")
  system(sprintf('cp %s %s/sorttable.js',st,outdir))
  system(sprintf('cp %s %s/dragtable.js',dt,outdir))
  for (j in names(mygs$index))
  {
    if ((nrow(mygs$res))>0)
    {
      intgenes <- intersect(as.character(detable$symbol),as.character(mygs$index[[j]]))
      ##idx1 <- unlist(lapply(strsplit(i,'_vs_'),function(x) x[1]))
      ##idx2 <- unlist(lapply(strsplit(i,'_vs_'),function(x) x[2]))
      xout <- detable[detable$symbol %in% intgenes,]
      ##idx1 <- colnames(xout)[grep(idx1,colnames(xout))] # Query
      ##idx2 <- colnames(xout)[grep(idx2,colnames(xout))] # Control
      normcols <- colnames(xout)[grep('normCounts',colnames(xout))]
      #                        selcols <- c('symbol','names','space','start','end','width','strand',normcols,'baseMean','log2FoldChange','stat','pvalue','padj','rej')
      xout <- xout[,selcols]
      xout <- xout[order(abs(xout[,vorder]),decreasing=TRUE),]
      isnum <- which(sapply(xout,is.numeric))
      for (coln in isnum) xout[,coln] <- round(xout[,coln],3)
      xout[is.na(xout)] <- 'NA'
      ##print(colnames(xout))
      imgsrc <- sprintf('../../../../../../%s/figs/stripcharts/%s.png',basename(DEdir),make.names(xout$names))
      ##xout$Plot <- sprintf('<a href="%s"><img src="%s" height=150 width=150></img></a>',imgsrc,imgsrc)
      ##xout$entrez <- sprintf('<a href="https://www.ncbi.nlm.nih.gov/gene/?term=%s">%s</a>',xout$entrez,xout$entrez)
      xout$symbol <- sprintf('<a href="http://www.informatics.jax.org/quicksearch/summary?queryType=exactPhrase&query=%s&submit=Quick+Search">%s</a>',xout$symbol,xout$symbol)
      fout <- file.path(outdir, sprintf('%s_genes.html',j))
      title <- sprintf('<center><h4>%s | %s: Genes in pathway %s\n Differential expression results</h4></center></p>',gs,i,j)
      p <- openPage(basename(fout),dirname(fout),title=title,link.javascript=c('sorttable.js','dragtable.js'),link.css='../../../../../../biostats.css')
      pp <- hwrite(title,i)
      pp <- paste(pp,hwrite(xout,table.class=list('sortable draggable'),table.align='center',row.names=FALSE),sep='\n')
      hwrite(pp,p)
      closePage(p)
    }
  }
}

#' Core HTML Assembler and Dynamic Table Generation Wrapper
#'
#' Writes dynamic table elements, formats embedded hyperlinks/images, and
#' outputs unified files utilizing internal `hwriter` constructs.
#'
#' @param x A structured data.frame to print.
#' @param file Character. Full path configuration target for file generation.
#' @param links A structured list containing column-matched URL targets.
#' @param tiny.pic A structured list containing column-matched graphic paths.
#' @param sorttable Character string containing contents of `sorttable.js`.
#' @param dragtable Character string containing contents of `dragtable.js`.
#' @param css Character. Relative path configuration target for CSS files.
#' @param title Character. Text element to serve as the document title.
#'
#' @export
write.html.mod2 <- function(x, file = paste0(htmlpath, htmlname), links = links,
                            tiny.pic = plots, sorttable = sorttable,
                            dragtable = dragtable,css='../../../../../biostats.css',title=title)
{
  ## Apply links and plots
  sel.links <- which(!sapply(links,is.null))
  sel.plots <- which(!sapply(tiny.pic,is.null))
  ans <- x
  for (j in sel.links) ans[,j] <- sprintf('<a href="%s">%s</a>',links[[j]],ans[,j])
  for (j in sel.plots) ans[,j] <- sprintf('<a href="%s"><img src="%s" height="100"></a>',tiny.pic[[j]],tiny.pic[[j]])
  ##ans
  ## Write .js files
  writeLines(sorttable,file.path(dirname(file),'sorttable.js'))
  writeLines(dragtable,file.path(dirname(file),'dragtable.js'))
  ## Write main output table
  p <- openPage(basename(file),dirname(file),link.javascript=c('sorttable.js','dragtable.js'),link.css='../../../../../biostats.css')
  pp <- hwrite(title)
  pp <- paste(pp,hwrite(ans,table.class=list('sortable draggable'),table.align='center',row.names=FALSE),sep='\n')
  hwrite(pp,p)
  closePage(p)
}

#' Low-Level ROAST Plot Engine and Property Structuring Adapter
#'
#' Generates functional PNG diagrams (GSEA, Stats, Heatmaps) from `roastgsa` data structures,
#' configures directory trees, and maps relational plotting variables to `write.html.mod2`.
#'
#' @param obj An object of class `roastgsa`.
#' @param htmlpath Character. Core output path where table indices will live.
#' @param htmlname Character. Main target sheet filename string.
#' @param plotpath Character. Relative path directory snippet where generated PNG files are exported.
#' @param plotstats Logical. If TRUE, renders distribution stats diagrams. Default TRUE.
#' @param plotgsea Logical. If TRUE, renders cumulative GSEA running sums. Default TRUE.
#' @param indheatmap Logical. If TRUE, renders expression matrix heatmaps. Default TRUE.
#' @param ploteffsize Logical. If TRUE, renders effect size metrics signature charts. Default TRUE.
#' @param links_plots A list containing explicit, override filepaths/URLs for generated plots.
#' @param y A numeric matrix of normalized gene expression values.
#' @param whplots Character vector. Explicit subset of pathway IDs to generate plots for.
#' @param geneDEhtmlfiles Character vector. Paths to matching child sheets created by `roastHtmlDETable`.
#' @param title Character. Text or HTML markup string to serve as header banner title.
#' @param margins Numeric vector of length 2. Axis margin controls for traditional heatmap rendering. Default c(15, 12).
#' @param sizesHeatmap Numeric vector of length 2. Canvas dimensions (Height, Width) for heatmap PNG exports. Default c(1200, 800).
#' @param typeheatmap Character vector. Target layout framework option string ("heatmap.2" vs "ggplot2").
#' @param intvar Character. Metadata factor tracking biological variable of interest.
#' @param adj.var Character. Design factor metadata variations to adjust for. Default NULL.
#' @param mycol Vector of color spectrum palettes for visual plotting elements.
#' @param varrot Numeric vector or matrix representing rotatable model variation coordinates.
#' @param psel Numeric vector. Extracted feature row pointer configurations. Default NULL.
#' @param sorttable Character string containing contents of `sorttable.js`.
#' @param dragtable Character string containing contents of `dragtable.js`.
#' @param ... Optional arguments passed down to downstream low-level plotting functions.
#'
#' @export
htmlrgsa2 <- function (obj, htmlpath = "", htmlname = "file.html", plotpath = "",
                       plotstats = TRUE, plotgsea = TRUE, indheatmap = TRUE, ploteffsize = TRUE,
                       links_plots = list(stats = NULL, gsea = NULL, heatmap = NULL,
                                          effsize = NULL), y, whplots = NULL, geneDEhtmlfiles = NULL,
                       title = "", margins = c(15, 12), sizesHeatmap = c(1200, 800),
                       typeheatmap = c("heatmap.2", "ggplot2"), intvar, adj.var = NULL,
                       mycol, varrot, psel = NULL, sorttable, dragtable, ...)
{
  if (!inherits(obj, "roastgsa"))
    stop("not a roastgsa object")
  if (ploteffsize)
    if (missing(varrot))
      stop("varrot is missing")
  x <- data.frame(geneset = rownames(obj$res), obj$res)
  index <- obj$index[rownames(x)]
  psel2 <- psel
  if (plotstats | plotgsea | indheatmap) {
    dir.create(paste0(htmlpath, plotpath))
    if (is.null(whplots))
      whplots <- names(index)
    if (!is.na(whplots[1])) {
      stats <- sort(obj$stats)
      index <- sapply(obj$index, function(x) which(names(stats) %in%
                                                     x))
      for (k in whplots) {
        if (plotstats) {
          png(paste0(htmlpath, plotpath, gsub("[[:punct:]]",
                                              " ", k), "_stats.png"))
          plotStats(obj, whplot = k, ...)
          dev.off()
        }
        if (plotgsea) {
          png(paste0(htmlpath, plotpath, gsub("[[:punct:]]",
                                              " ", k), "_gsea.png"))
          plotGSEA(obj, whplot = k, ...)
          dev.off()
        }
        if (indheatmap) {
          png(paste0(htmlpath, plotpath, gsub("[[:punct:]]",
                                              " ", k), "_heatmap.png"), width = sizesHeatmap[2],
              height = sizesHeatmap[1])
          if (typeheatmap[1] == "ggplot2")
            heatmaprgsa_hm(obj, y, whplot = k, mycol = mycol,
                           intvar = intvar, adj.var = adj.var, psel = psel2,
                           ...)
          else heatmaprgsa_hm(obj, y, whplot = k, mycol = mycol,
                              intvar = intvar, adj.var = adj.var, psel = psel2,
                              ...)
          dev.off()
        }
        if (ploteffsize) {
          png(paste0(htmlpath, plotpath, gsub("[[:punct:]]",
                                              " ", k), "_effsize.png"))
          ploteffsignaturesize(obj, varrot, whplot = k)
          dev.off()
        }
      }
    }
  }
  if (!is.null(geneDEhtmlfiles))
    x$geneDEinfo <- rep("view", dim(x)[1])
  if (plotstats)
    x$plot_stats <- NA
  if (plotgsea)
    x$plot_gsea <- NA
  if (indheatmap)
    x$heatmap <- NA
  if (ploteffsize)
    x$plot_effsize <- NA
  links <- vector("list", length = ncol(x))
  names(links) <- colnames(x)
  plots <- links
  if (!is.null(links_plots$stats))
    links$plot_stats <- plots$plot_stats <- links_plots$stats
  else {
    if (plotstats)
      links$plot_stats <- plots$plot_stats <- paste0(plotpath,
                                                     gsub("[[:punct:]]", " ", rownames(x)), "_stats.png")
  }
  if (!is.null(links_plots$gsea))
    links$plot_gsea <- plots$plot_gsea <- links_plots$gsea
  else {
    if (plotgsea)
      links$plot_gsea <- plots$plot_gsea <- paste0(plotpath,
                                                   gsub("[[:punct:]]", " ", rownames(x)), "_gsea.png")
  }
  if (!is.null(links_plots$heatmap))
    links$heatmap <- plots$heatmap <- links_plots$heatmap
  else {
    if (indheatmap)
      links$heatmap <- plots$heatmap <- paste0(plotpath,
                                               gsub("[[:punct:]]", " ", rownames(x)), "_heatmap.png")
  }
  if (!is.null(links_plots$heatmap))
    links$effsize <- plots$effsize <- links_plots$effsize
  else {
    if (ploteffsize)
      links$plot_effsize <- plots$plot_effsize <- paste0(plotpath,
                                                         gsub("[[:punct:]]", " ", rownames(x)), "_effsize.png")
  }
  if (!is.null(geneDEhtmlfiles)) {
    links$geneDEinfo <- rep(NA, dim(x)[1])
    links$geneDEinfo[1:length(geneDEhtmlfiles)] <- geneDEhtmlfiles
  }
  write.html.mod2(x, file = paste0(htmlpath, htmlname), links = links,
                  tiny.pic = plots, title = title, sorttable = sorttable,
                  dragtable = dragtable, ...)
}

## The End...
