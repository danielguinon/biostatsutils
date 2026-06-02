#' Generate Linear Contrast Vectors from Model Formula and Group Selectors
#'
#' This function evaluates group subsets from a data frame based on user-defined
#' criteria, computes the design matrix means for those groups (either balanced
#' across factor levels or based on observed data), and calculates the specified
#' linear contrast.
#'
#' @param d A data frame containing the design variables and covariates.
#' @param form A formula object (e.g., \code{y ~ x1 + x2}) specifying the model design.
#' @param lsel A named list of lists defining the groups. Each sub-list contains
#'   variable-value pairs used to subset the data (e.g., \code{list(Group1 = list(sex = "M"))}).
#' @param signs A numeric vector of weights/signs indicating how to contrast the
#'   groups (e.g., \code{c(1, -1)}).
#' @param type A character string specifying the weighting type; either
#'   \code{'balanced'} (default) to average over grid levels, or another value to
#'   account for observed means.
#' @param ngroups An optional character vector providing custom names for the defined groups.
#'
#' @return A list containing three elements:
#' \block{
#'   \item{grs}{A matrix of the average design matrix rows for each defined group.}
#'   \item{cont}{The computed linear contrast vector (calculated via matrix multiplication of \code{grs} and \code{signs}).}
#'   \item{sels}{A logical matrix indicating which rows of the original data frame \code{d} belong to each group selector.}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' # Example setup (assuming a data frame 'df' with factor 'treatment')
#' my_formula <- response ~ treatment
#' group_selectors <- list(T1 = list(treatment = "A"), T2 = list(treatment = "B"))
#' contrast_weights <- c(1, -1)
#'
#' make.contrasts(d = df, form = my_formula, lsel = group_selectors, signs = contrast_weights)
#' }
make.contrasts <- function(d, form, lsel, signs, type='balanced', ngroups=NULL)
{

  nv <- strsplit(paste(form)[2], split=' \\+ ')[[1]];
  nv <- nv[regexpr("\\:", nv) < 0];
  nv <- nv[nv!='-1']
  for (n in nv) if (is.character(d[, n])) d[, n] <- factor(d[, n])

  nvf <- nv[sapply(nv, function(v, d) is.factor(d[, v]), d)];
  nvn <- nv[!nv%in%nvf];

  daf <- expand.grid(sapply(nvf, function(o) levels(d[, o]), simplify=F));
  dan <- matrix(rep(apply(d[, nvn, drop=F], 2, mean), each=nrow(daf)), nrow=nrow(daf))
  colnames(dan) <- nvn;
  colnames(daf) <- nvf;
  dal <- cbind(daf, dan)[nv];

  dsg <- model.matrix(form, dal);

  lss <- sapply(lsel, function(sel, dal)
  {
    ss <- sapply(1:length(sel), function(j, sel, dal)
    {
      nm <- names(sel)[[j]];
      if (is.factor(dal[, nm])) paste("'", sel[[j]], "'", sep="");
    }, sel, dal, simplify=F);
    names(ss) <- names(sel);
    ss;
  }, dal, simplify=F);

  ds0 <- sapply(lss, function(ss, d)
  {
    eval(parse(text=paste(sapply(1:length(ss), function(j, ss, d)
    {
      o <- ss[[j]];
      names(o) <- names(ss)[j];
      paste("d$", names(ss)[j], "%in%c(",
            paste(o, collapse=', ', sep=''), ")", sep='');
    },  ss, dal), collapse=' & ')))
  }, d)

  if (any(apply(ds0, 1, sum) > 1)) stop("Error: Overlapping groups!");

  if (type=='balanced')
  {
    ds <- sapply(lss, function(ss, dal)
    {
      eval(parse(text=paste(sapply(1:length(ss), function(j, ss, dal)
      {
        o <- ss[[j]];
        names(o) <- names(ss)[j];
        paste("dal$", names(ss)[j], "%in%c(",
              paste(o, collapse=', ', sep=''), ")", sep='');
      },  ss, dal), collapse=' & ')))
    }, dal)

    grs <- t(apply(ds, 2, function(s, dsg) apply(dsg[s, , drop=F], 2, mean), dsg));
    cont <- (t(grs)%*%signs)[, 1]

  }else
  {
    dsg0 <- model.matrix(form, d);

    nvm <- colnames(dsg)[regexpr("\\:", colnames(dsg)) < 0];
    lnve <- sapply(ldsg, function(ds, nvm)
    {
      nvm[apply(ds[, nvm, drop=F], 2, function(o) length(unique(o)) == 1)];
    }, nvm);
    lnvd <- sapply(ldsg, function(ds, nvm)
    {
      nvm[apply(ds[, nvm, drop=F], 2, function(o) length(unique(o)) > 1)];
    }, nvm);
    lmed.nvd <- sapply(lnvd, function(v, dsg0) apply(dsg0[, v, drop=F], 2, mean), dsg0);

    lp <- sapply(1:length(med.nvd), function(j, nvm)
    {
      p <- data.frame(t(c(unique(ldsg[[j]][, lnve[[j]]]), lmed.nvd[[j]])));
      colnames(p) <- c(lnve[[j]], lnvd[[j]]);
      p[, nvm, drop=F];
    }, nvm);

    grs <- t(apply(ds, 2, function(s, dsg) apply(dsg[s, , drop=F], 2, mean), dsg));
    cont <- (t(grs)%*%signs)[, 1]

  }

  if (!is.null(ngroups)) rownames(grs) <- colnames(ds0) <- ngroups;
  list(grs=grs, cont=cont, sels=ds0);

}
