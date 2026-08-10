#' Execute ROAST Gene Set Analysis and Output Nested HTML Reports
#'
#' Evaluates multiple collections of gene sets across multiple contrasts using
#' parallelized ROAST operations, generating structured downstream HTML summaries.
#'
#' @param x A list of lists containing \code{roastgsa} results objects structured by contrast and gene set collection.
#' @param y A list of data.frames containing Differential Expression (DE) metrics structured by contrast.
#' @param mat A numeric matrix of normalized gene expression values.
#' @param mc.cores Default 1.
#' @param outdir Character. Output path for results directories. Default './'.
#' @param maxgs Integer. Maximum number of top paths to include per report. Default 50.
#' @param indhtml Logical. If TRUE, creates dedicated nested child HTML sheets for individual genes. Default TRUE.
#' @param DEdir Character. Source directory path containing static figures (e.g., stripcharts). Default NULL.
#' @param returnData Logical. If TRUE, returns the nested results list structure. Default TRUE.
#' @param intvar Character. Name of the primary design variable of interest. Default NULL.
#' @param mycol Vector of color definitions for heatmap generation.
#'
#' @importFrom parallel mclapply
#' @export
roastHtmlTables <- function (x, y, mat, mc.cores = 1, outdir = "./",
  maxgs = 50, indhtml = TRUE, DEdir = NULL,
  returnData = TRUE, intvar = NULL, mycol = NULL){
  if (is.null(mycol))
    stop("A color palette 'mycol' must be provided.")
  if (!dir.exists(outdir))
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  sorttable_src <- system.file("javascript", "sorttable.js", package = "phenoTest")
  dragtable_src <- system.file("javascript", "dragtable.js", package = "phenoTest")

  if (sorttable_src == "" || dragtable_src == "") {
    warning("JavaScript utility assets could not be located via system.file. Falling back to empty strings.")
    sorttable_code <- ""
    dragtable_code <- ""
  } else {
    sorttable_code <- readLines(sorttable_src, warn = FALSE)
    dragtable_code <- readLines(dragtable_src, warn = FALSE)
  }

  # All combinations
  grid_tasks <- expand.grid(
    gs = names(x[[1]]),
    i  = names(y),
    stringsAsFactors = FALSE
  )

  results_flat <- mclapply(seq_len(nrow(grid_tasks)), function(idx) {
    gs <- grid_tasks$gs[idx]
    i  <- grid_tasks$i[idx]

    mygs <- x[[i]][[gs]]
    geneset_dir <- file.path(outdir, "roastGSA", "html", gs)
    filename_html <- sprintf("roastGSA_MaxMean_%s_%s.html", gs, i)

    roastHtmlTable(
      mygs = mygs, gs = gs, i = i, filename = filename_html,
      out_dirname = geneset_dir, indhtml = indhtml, DEdir = DEdir,
      detable = y[[i]], maxgs = maxgs, mat = mat, intvar = intvar,
      selcols = colnames(y[[i]]), vorder = colnames(y[[i]])[ncol(y[[i]])],
      mycol = mycol, outdir = outdir, sorttable = sorttable_code,
      dragtable = dragtable_code
    )
  }, mc.cores = mc.cores)

  if (returnData) {
    ans <- split(results_flat, grid_tasks$gs)
    ans <- lapply(ans, function(gs_group) {
      names(gs_group) <- grid_tasks$i[1:length(gs_group)]
      gs_group
    })
    return(ans)
  }
}


#' Generate HTML Report Table for a Specific Gene Set and Contrast
#'
#' Filters gene set analysis results based on significance controls, builds underlying
#' image assets, and generates a structured summary web page directory tracking background signals.
#'
#' @param mygs A specific \code{roastgsa} results object.
#' @param gs Character. The name or identifier of the gene set collection.
#' @param i Character. The name or identifier of the contrast group.
#' @param filename Character. Target file name string for the output layout page.
#' @param dirname Character. Directory path destination tracking framework components.
#' @param indhtml Logical. If TRUE, generates nested child pages for local gene expressions.
#' @param DEdir Character. Optional source path containing auxiliary diagnostic illustrations. Default NULL.
#' @param detable A data.frame containing differential expression parameters for subsetting workflows. Default NULL.
#' @param maxgs Integer. Maximum threshold cutoff for tracking high-scoring gene expressions. Default 50.
#' @param verbose Logical. If TRUE, messages execution metrics out to console. Default TRUE.
#' @param mat A numeric matrix of expression value counts.
#' @param intvar Character. Primary experimental variable context tracking group metadata. Default NULL.
#' @param selcols Vector of string characters matching explicit column subsets to extract.
#' @param vorder Character. Destination column sorting sequence assignment identifier.
#' @param apval.cut Numeric. Adjusted p-value threshold for considering paths significant. Default 0.05.
#' @param mycol Vector of color definitions for standard mapping layouts. Default NULL.
#' @param outdir Character. Base directory tracking root structural framework components. Default './'.
#' @param sorttable Character string containing raw JavaScript file code logic for sortable structures. Default "".
#' @param dragtable Character string containing raw JavaScript file code logic for draggable tables. Default "".
#'
#' @export
roastHtmlTable <- function(mygs, gs, i, filename, out_dirname, indhtml, DEdir = NULL,
                           detable = NULL, maxgs = 50, verbose = TRUE, mat, intvar = NULL,
                           selcols, vorder, apval.cut = 0.05, mycol = NULL, outdir = "./",
                           sorttable = "", dragtable = "")
{
  if (indhtml && is.null(detable))
    stop("indhtml is TRUE but no DE table provided")

  # Guard: ensure result frame exists and is non-empty
  if (is.null(mygs$res) || nrow(mygs$res) == 0) {
    if (verbose) message(sprintf("Skipping %s in %s: empty results.", gs, i))
    return(NULL)
  }

  # Sort by absolute NES decreasing
  mygs$res <- mygs$res[order(abs(mygs$res$nes), decreasing = TRUE), ]

  # select significant gene sets (handling NAs)
  sel <- which(mygs$res$adj.pval < apval.cut)
  if (length(sel) < 10) {
    mmax <- min(10, nrow(mygs$res))
    sel <- seq_len(mmax)
  }

  mygs$res <- mygs$res[sel, , drop = FALSE]

  # Trim to maxgs
  if (nrow(mygs$res) > 0) {
    mygs$res <- mygs$res[seq_len(min(nrow(mygs$res), maxgs)), , drop = FALSE]
  }

  # Align gene set indices
  mygs$index <- mygs$index[rownames(mygs$res)]

  if (verbose) {
    message(sprintf("Writing output for %s in %s (%d selected genesets)",
                    gs, i, nrow(mygs$res)))
  }

  if (!dir.exists(out_dirname))
    dir.create(out_dirname, recursive = TRUE, showWarnings = FALSE)


  if (indhtml) {
    geneDEhtmlfiles <- sprintf("%s_indhtml/%s_genes.html", i, rownames(mygs$res))
  } else {
    geneDEhtmlfiles <- NULL
  }

  html_target_path <- file.path(outdir, "roastGSA", "html", gs)

  htmlrgsa2(
    obj = mygs,
    htmlpath = html_target_path,
    htmlname = filename,
    plotpath = sprintf("%s_images/", i),
    indheatmap = TRUE,
    y = mat,
    intvar = intvar,
    ploteffsize = FALSE,
    mycol = mycol,
    geneDEhtmlfiles = geneDEhtmlfiles,
    sorttable = sorttable,
    dragtable = dragtable,
    whplot = rownames(mygs$res),
    title = sprintf("<center><h4>%s | %s</h4></center>", gs, i)
  )

  if (indhtml) {
    roastHtmlDETable(
      mygs = mygs, gs = gs, i = i, detable = detable,
      DEdir = DEdir, filename = filename, out_dirname = out_dirname,
      maxgs = maxgs, verbose = verbose, selcols = selcols,
      vorder = vorder
    )
  }
}


#' Export Nested Differential Expression Target Layout Sub-tables
#'
#' Evaluates specific matching features between target gene set lists and dynamic sequence structures,
#' writing localized independent HTML tracking sheets mapping data frames into interactive sortable models.
#'
#' @param mygs A specific \code{roastgsa} results structure object.
#' @param gs Character. The target name context matching background collections.
#' @param i Character. The specific contrast assignment tracker.
#' @param detable A data.frame mapping baseline differential parameters.
#' @param DEdir Character. Input file storage directory referencing static plots.
#' @param filename Character. Output name specification tracking target templates.
#' @param dirname Character. Container folder path destination mapping components.
#' @param maxgs Integer. Maximum scale tracking threshold. Default 50.
#' @param verbose Logical. Enables output standard messages. Default TRUE.
#' @param selcols Vector of characters filtering acceptable column outputs.
#' @param vorder Character string declaring variable order hierarchy sequence targeting.
#'
#' @importFrom hwriter hwrite openPage closePage
#' @export
roastHtmlDETable <- function(mygs, gs, i, detable, DEdir, filename, out_dirname,
                             maxgs=50, verbose=TRUE, selcols, vorder) {

  # if (!"symbol" %in% colnames(detable)) {
  #   stop("The provided differential expression table 'detable' is missing a 'symbol' column.")
  # }
  if (!"symbol" %in% colnames(detable)) {
    warning(sprintf("Skipping gene-level DE tables for %s | %s: no 'symbol' column (found: %s).",
                    gs, i, paste(colnames(detable), collapse = ", ")))
    return(invisible(NULL))

  }

  outdir <- file.path(out_dirname, sprintf('%s_indhtml', i))
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  st <- system.file("javascript", "sorttable.js", package = "phenoTest")
  dt <- system.file("javascript", "dragtable.js", package = "phenoTest")

  if (st != "") file.copy(st, file.path(outdir, 'sorttable.js'), overwrite = TRUE)
  if (dt != "") file.copy(dt, file.path(outdir, 'dragtable.js'), overwrite = TRUE)

  for (j in names(mygs$index)) {
    if ((nrow(mygs$res)) > 0) {
      intgenes <- intersect(as.character(detable$symbol), as.character(mygs$index[[j]]))
      xout <- detable[detable$symbol %in% intgenes, ]

      existing_cols <- intersect(selcols, colnames(xout))
      xout <- xout[, existing_cols, drop = FALSE]

      if (vorder %in% colnames(xout)) {
        xout <- xout[order(abs(xout[, vorder]), decreasing = TRUE), ]
      }

      isnum <- which(sapply(xout, is.numeric))
      for (coln in isnum) xout[, coln] <- round(xout[, coln], 3)
      xout[is.na(xout)] <- 'NA'

      xout$symbol <- sprintf('<a href="http://www.informatics.jax.org/quicksearch/summary?queryType=exactPhrase&query=%s&submit=Quick+Search">%s</a>', xout$symbol, xout$symbol)

      fout <- file.path(outdir, sprintf('%s_genes.html', j))
      title <- sprintf('<center><h4>%s | %s: Genes in pathway %s\n Differential expression results</h4></center></p>', gs, i, j)

      p <- openPage(basename(fout), dirname(fout), title = title,
                    link.javascript = c('sorttable.js', 'dragtable.js'),
                    link.css = '../../../../../../biostats.css')
      pp <- hwrite(title, i)
      pp <- paste(pp, hwrite(xout, table.class = list('sortable draggable'), table.align = 'center', row.names = FALSE), sep = '\n')
      hwrite(pp, p)
      closePage(p)
    }
  }
}


#' Generate Global Enrichment Tables and Structural Diagnostic Layout Graphics
#'
#' Processes statistical records across designated parameters to generate structural png
#' matrix diagnostics and wraps summary matrix metadata maps into primary target pages.
#'
#' @param obj A valid \code{roastgsa} container metrics object.
#' @param htmlpath Character. Core framework folder directory target path. Default "".
#' @param htmlname Character. Main target document layout filename. Default "file.html".
#' @param plotpath Character. Local child relative folder directory targeting internal visuals. Default "".
#' @param plotstats Logical. Controls running diagnostic visual statistics generation loops. Default TRUE.
#' @param plotgsea Logical. Controls execution pathways mapping baseline GSEA profiles. Default TRUE.
#' @param indheatmap Logical. Renders distinct matrix heat maps tracking active target expressions. Default TRUE.
#' @param ploteffsize Logical. Evaluates profile size effects signatures. Default TRUE.
#' @param links_plots A structured list parsing unique relative paths matching downstream plot destinations.
#' @param y Expression counts target numerical background matrix.
#' @param whplot Target array matching identifiers across gene set indexes to render manually. Default NULL.
#' @param geneDEhtmlfiles Character mapping array pointing cleanly to lower-level nested child page frames. Default NULL.
#' @param title Character layout header text string mapping raw web content templates. Default "".
#' @param margins Numeric vector tracking visual boundaries formatting parameters. Default c(15, 12).
#' @param sizesHeatmap Numeric vector parsing layout dimension specifications. Default c(1200, 800).
#' @param typeheatmap Character vector detailing preferred baseline engine types. Default c("heatmap.2", "ggplot2").
#' @param intvar Character parameter context specifying active experimental design tracking parameters.
#' @param adj.var Supplementary metadata tracking adjustments covariates matrix. Default NULL.
#' @param mycol Vector palette configuring structural map visualization themes.
#' @param varrot Numeric signature parameter targeting profile metrics calculations. Optional.
#' @param psel Supplementary significance parameter maps selector. Default NULL.
#' @param sorttable String containing base Javascript utility logic. Default "".
#' @param dragtable String containing baseline Javascript interactivity logic. Default "".
#' @param ... Additional argument options passed securely to downstream graphics functions.
#'
#' @importFrom grDevices png dev.off
#' @export
htmlrgsa2 <- function (obj, htmlpath = "", htmlname = "file.html", plotpath = "",
                       plotstats = TRUE, plotgsea = TRUE, indheatmap = TRUE, ploteffsize = TRUE,
                       links_plots = list(stats = NULL, gsea = NULL, heatmap = NULL, effsize = NULL),
                       y, whplot = NULL, geneDEhtmlfiles = NULL, title = "", margins = c(15, 12),
                       sizesHeatmap = c(1200, 800), typeheatmap = c("heatmap.2", "ggplot2"),
                       intvar, adj.var = NULL, mycol, varrot, psel = NULL,
                       sorttable = "", dragtable = "", ...)
{
  if (!inherits(obj, "roastgsa"))
    stop("not a roastgsa object")
  if (ploteffsize && missing(varrot))
    stop("varrot is missing and required when ploteffsize is TRUE")

  htmlpath <- sub("/*$", "/", htmlpath) # This line guarantees one trailing slash "/"

  # Unique, filesystem-safe id per gene set — reused everywhere a filename
  # is derived from a gene-set name, so collisions can't silently overwrite plots.
  x <- data.frame(geneset = rownames(obj$res), obj$res)
  safe_id <- setNames(
    sprintf("%04d_%s", seq_along(rownames(x)), gsub("[^[:alnum:]]+", "_", rownames(x))),
    rownames(x)
  )

  index <- obj$index[rownames(x)]
  psel2 <- psel

  if (plotstats | plotgsea | indheatmap) {
    dir.create(paste0(htmlpath, plotpath), recursive = TRUE, showWarnings = FALSE)
    if (is.null(whplot))
      whplot <- names(index)

    if (!is.na(whplot[1])) {
      stats <- sort(obj$stats)
      # index <- sapply(obj$index, function(z) which(names(stats) %in% z))

      for (k in whplot) {
        # clean_k <- gsub("[[:punct:]]", " ", k)
        clean_k <- safe_id[[k]]

        # if (plotstats) {
        #   png(paste0(htmlpath, plotpath, clean_k, "_stats.png"))
        #   plotStats(obj, whplot = k, ...)
        #   dev.off()
        # }
        if (plotstats) {
          png(paste0(htmlpath, plotpath, clean_k, "_stats.png"))
          tryCatch(plotStats(obj, whplot = k, ...),
                   error = function(e) message("plotStats failed for '", k, "': ", conditionMessage(e)))
          dev.off()
        }
        if (plotgsea) {
          png(paste0(htmlpath, plotpath, clean_k, "_gsea.png"))
          tryCatch(plotGSEA(obj, whplot = k, ...),
                   error = function(e) message("plotStats failed for '", k, "': ", conditionMessage(e)))
          dev.off()
        }
        if (indheatmap) {
          png(paste0(htmlpath, plotpath, clean_k, "_heatmap.png"),
              width = sizesHeatmap[2], height = sizesHeatmap[1])
          tryCatch({
            hm_obj <- heatmaprgsa_hm(obj, y = y, whplot = k, mycol = mycol,
                                     intvar = intvar, adj.var = adj.var, psel = psel2,
                                     toplot = TRUE, pathwaylevel = FALSE, ...)
            if(!is.null(hm_obj) && inherits(hm_obj, c("ggplot", "gg", "grob", "trellis"))) print(hm_obj)
          }, error = function(e) message("heatmaprgsa_hm failed for '", k, "': ", conditionMessage(e)))
          dev.off()
        }
        if (ploteffsize) {
          png(paste0(htmlpath, plotpath, clean_k, "_effsize.png"))
          tryCatch(ploteffsignaturesize(obj, whplot = k, ...),
                   error = function(e) message("plotStats failed for '", k, "': ", conditionMessage(e)))
          dev.off()
        }
      }
    }
  }

  if (!is.null(geneDEhtmlfiles))
    x$geneDEinfo <- rep("view", dim(x)[1])
  if (plotstats)  x$plot_stats <- NA
  if (plotgsea)   x$plot_gsea <- NA
  if (indheatmap) x$heatmap <- NA
  if (ploteffsize) x$plot_effsize <- NA

  links <- vector("list", length = ncol(x))
  names(links) <- colnames(x)
  plots <- links

  # clean_rownames <- gsub("[[:punct:]]", " ", rownames(x))
  clean_rownames <- safe_id[rownames(x)]


  if (!is.null(links_plots$stats)) links$plot_stats <- plots$plot_stats <- links_plots$stats
  else if (plotstats) links$plot_stats <- plots$plot_stats <- paste0(plotpath, clean_rownames, "_stats.png")

  if (!is.null(links_plots$gsea)) links$plot_gsea <- plots$plot_gsea <- links_plots$gsea
  else if (plotgsea) links$plot_gsea <- plots$plot_gsea <- paste0(plotpath, clean_rownames, "_gsea.png")

  if (!is.null(links_plots$heatmap)) links$heatmap <- plots$heatmap <- links_plots$heatmap
  else if (indheatmap) links$heatmap <- plots$heatmap <- paste0(plotpath, clean_rownames, "_heatmap.png")

  if (!is.null(links_plots$effsize)) links$effsize <- plots$effsize <- links_plots$effsize
  else if (ploteffsize) links$plot_effsize <- plots$plot_effsize <- paste0(plotpath, clean_rownames, "_effsize.png")

  if (!is.null(geneDEhtmlfiles)) {
    links$geneDEinfo <- rep(NA, dim(x)[1])
    links$geneDEinfo[1:length(geneDEhtmlfiles)] <- geneDEhtmlfiles
  }

  write.html.mod2(x, file = paste0(htmlpath, htmlname), links = links,
                  tiny.pic = plots, title = title, sorttable = sorttable,
                  dragtable = dragtable, ...)
}



#' Write Structured Web Documentation Files onto Persistent Target Disk Spaces
#'
#' Merges string lists mapping anchor tags and thumbnail targets directly into data frames,
#' writing target configurations alongside local dependencies onto persistent file assets.
#'
#' @param x A structured data.frame containing metadata layout results tracking values.
#' @param file Character. Clean absolute target file path to compile toward. Default "file.html".
#' @param links List parsing functional string formatting directives targeting anchor parameters. Default list().
#' @param tiny.pic List mapping target layout parameters defining visual image destinations. Default list().
#' @param sorttable String sequence detailing sorting utilities logic scripts. Default "".
#' @param dragtable String sequence tracking dragging utilities engine parameters. Default "".
#' @param css Character string declaring stylesheet layout configurations. Default '../../../../../biostats.css'.
#' @param title Character string header title injected onto page setups. Default "".
#'
#' @importFrom hwriter hwrite openPage closePage
#' @export
write.html.mod2 <- function(x, file = "file.html", links = list(),
                            tiny.pic = list(), sorttable = "",
                            dragtable = "", css = '../../../../../biostats.css', title = "")
{
  sel.links <- which(!sapply(links, is.null))
  sel.plots <- which(!sapply(tiny.pic, is.null))
  ans <- x

  for (j in sel.links) ans[,j] <- sprintf('<a href="%s">%s</a>', links[[j]], ans[,j])
  for (j in sel.plots) ans[,j] <- sprintf('<a href="%s"><img src="%s" height="100"></a>', tiny.pic[[j]], tiny.pic[[j]])

  if (!dir.exists(dirname(file))) dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)

  writeLines(sorttable, file.path(dirname(file), 'sorttable.js'))
  writeLines(dragtable, file.path(dirname(file), 'dragtable.js'))

  p <- openPage(basename(file), dirname(file), link.javascript = c('sorttable.js', 'dragtable.js'), link.css = css)
  pp <- hwrite(title)
  pp <- paste(pp, hwrite(ans, table.class = list('sortable draggable'), table.align = 'center', row.names = FALSE), sep = '\n')
  hwrite(pp, p)
  closePage(p)
}
