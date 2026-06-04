#' Write a Data Frame to an Interactive HTML Table
#'
#' @param x A data.frame.
#' @param links A list of URL character vectors matching `x` columns, or NULL.
#' @param tiny.pic A list of thumbnail image paths matching `x` columns, or NULL.
#' @param tiny.pic.size Numeric. Width and height of thumbnails in pixels. Default 100.
#' @param title Character. Table title.
#' @param file Character. Output file path.
#' @param digits Integer. Decimal places for numeric columns. Default 3.
#'
#' @export
write.html.acm <-
  function (x, links, tiny.pic, tiny.pic.size = 100, title = "",
            file, digits = 3)
  {
    stopifnot(class(x) == "data.frame")
    if (missing(links))
      links <- vector("list", ncol(x))
    if (missing(tiny.pic))
      tiny.pic <- vector("list", ncol(x))
    stopifnot(class(links) == "list")
    stopifnot(class(tiny.pic) == "list")
    stopifnot(length(links) == ncol(x))
    stopifnot(length(tiny.pic) == ncol(x))
    stopifnot(!missing(file))
    column.class <- unlist(lapply(x, class))
    for (j in 1:ncol(x)) {
      if (column.class[j] == "factor")
        x[, j] <- as.character(x[, j])
      if (column.class[j] == "numeric")
        x[, j] <- round(x[, j], digits = digits)
    }
    cat("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">\n",
        sep = "", file = file, append = F)
    cat("<html>\n", file = file, append = T)
    cat("<body>\n", file = file, append = T)
    cat(paste("<CAPTION ALIGN=\"top\"><center><B>", title, "</B></center></CAPTION><BR>\n"),
        sep = "", file = file, append = T)
    cat("<TABLE border=1>\n", file = file, append = T)
    cat("<TR>\n", file = file, append = T)
    for (j in 0:ncol(x)) {
      cat("<TH>", file = file, append = T)
      if(j == 0)
        cat("", file = file, append = T)
      else
        cat(colnames(x)[j], file = file, append = T)
      cat("</TH>\n", file = file, append = T)
    }
    cat("</TR>\n", file = file, append = T)
    for (i in 1:nrow(x)) {
      cat("<TR>\n", file = file, append = T)
      cat("<TH>", file = file, append = T)
      cat(rownames(x)[i], file = file, append = T)
      cat("</TH>\n", file = file, append = T)
      for (j in 1:ncol(x)) {
        cat("<TD>", file = file, append = T)
        if (is.null(links[[j]]) & is.null(tiny.pic[[j]])) {
          cat(x[i, j], file = file, append = T)
        }
        else if (is.null(links[[j]]) & !is.null(tiny.pic[[j]])) {
          cat(paste("<A HREF=\"", links[[j]][[i]], "\"><img src=\"",
                    tiny.pic[[j]][[i]], "\" height=\"", tiny.pic.size,
                    "\" width=\"", tiny.pic.size, "\" /></A>",
                    sep = ""), file = file, append = T)
        }
        else if (!is.null(links[[j]]) & is.null(tiny.pic[[j]])) {
          cat(paste("<A HREF=\"", links[[j]][[i]], "\">",
                    x[i, j], "</A>", sep = ""), file = file, append = T)
        }
        else if (!is.null(links[[j]]) & !is.null(tiny.pic[[j]])) {
          cat(paste("<A HREF=\"", links[[j]][[i]], "\"><img src=\"",
                    tiny.pic[[j]][[i]], "\" height=\"", tiny.pic.size,
                    "\" width=\"", tiny.pic.size, "\" /></A>",
                    sep = ""), file = file, append = T)
        }
        cat("</TD>\n", file = file, append = T)
      }
      cat("</TR>\n", file = file, append = T)
    }
    cat("</TABLE>\n", file = file, append = T)
    cat("</body>\n", file = file, append = T)
    cat("</html>\n", file = file, append = T)
    phenoTest:::sortDragHtmlTable(filename = file)
  }
