# Copyright 2019, 2020 Steve Yadlowsky
# This includes modifications to software provided under the following copyright.
#
#  Copyright (C) 1995-2015 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

fast_polym <- function (..., degree = 1, coefs = NULL, vals_per_dim = NULL, raw = FALSE)
{
  if (is.matrix(...)) {
    dots <- unclass(as.data.frame(...))
  }
  else {
    dots <- list(...)
  }
    nd <- length(if(is.null(coefs)) dots else coefs)  # number of variables
    if (nd == 0)
        stop("must supply one or more vectors")
    if (degree > nd) {
			## z:= combinations of (0:degree) of all variables
			z <- do.call(expand.grid,
									 c(rep.int(list(0:degree), nd), KEEP.OUT.ATTRS = FALSE))
			## sum of all degrees must be in  1:degree :
			s <- rowSums(z)
			ind <- 0 < s  &  s <= degree
			z <- z[ind, , drop=FALSE]
			s <- s[ind]
    } else {
      locs <- expand.grid(c(rep.int(list(0:nd), degree)))[-1,]
      if (degree > 1) {
        nondecreasing_locs <- apply(locs, 1, function(row) { all(row==cummax(row)) })
        locs <- locs[nondecreasing_locs,]
      } else {
        locs <- matrix(locs, ncol=1)
      }
      z <- matrix(0, ncol=nrow(locs), nrow=nd)
      for (i in 1:degree) {
        assignment <- locs[,i]
        degree_locs <- vapply(assignment, function(loc, len) { s<-rep(0, len); s[loc]<-1; s }, rep(0,nd), nd)
				z <- z + degree_locs
      }
      z <- t(z)
      s <- rowSums(z)
    }
    if (is.null(vals_per_dim)) {
      vals_per_dim <- vapply(dots, function(feature) { length(unique(feature)) }, 1)
    }
    more_vals_than_degree <- apply(z, 1, function(row) {all(row < vals_per_dim)})
    z <- z[more_vals_than_degree,,drop=F]
    if(is.null(coefs)) {
	aPoly <- poly(dots[[1L]], pmin(degree, vals_per_dim[1]-1),
		      raw = raw, simple = raw && nd > 1)
	if (nd == 1)
	    return(aPoly)
	## nd >= 2 from here on
	n <- lengths(dots)
	if (any(n != n[1L]))
	    stop("arguments must have the same length")
	res <- cbind(1, aPoly)[, 1L +z[, 1]]
	## attribute "coefs" = list of coefs from individual variables
	if (!raw) coefs <- list(attr(aPoly, "coefs"))
	for (i in 2:nd) {
	    aPoly <- poly(dots[[i]], pmin(degree, vals_per_dim[i]-1),
			  raw = raw, simple = raw)
	    res <- res * cbind(1, aPoly)[, 1L +z[, i]]
	    if (!raw) coefs <- c(coefs, list(attr(aPoly, "coefs")))
	}
	colnames(res) <- apply(z, 1L, function(x) paste(x, collapse = "."))
	structure(res,
		  degree =  as.vector(s),
		  vals_per_dim = vals_per_dim,
		  coefs = if (!raw) coefs,
		  class = c("fastpoly", "poly", "matrix"))
    } else { ## use 'coefs' for prediction
	newdata <- as.data.frame(dots) # new data
	if (nd != ncol(newdata))
	    stop("wrong number of columns in new data: ", deparse(substitute(...)))
	res <- cbind(1, poly(newdata[[1]], degree=pmin(degree, vals_per_dim[1]-1),
			     coefs=coefs[[1]], simple=TRUE))[, 1L +z[, 1]]
	if(nd > 1) for (i in 2:nd)
            res <- res*cbind(1, poly(newdata[[i]], degree=pmin(degree, vals_per_dim[i]-1),
                                     coefs=coefs[[i]], simple=TRUE))[, 1L +z[, i]]
	colnames(res) <- apply(z, 1L, function(x) paste(x, collapse = "."))
        ## no 'coefs' and 'degree', nor "poly" class
	res
    }
}


predict.fastpoly <- function(object, newdata, ...)
{
    if(missing(newdata)) {
	object
    }
    else if(is.null(attr(object, "coefs"))) {
	fast_polym(newdata, degree = max(attr(object, "degree")),
             raw = TRUE)
    }
    else {
	fast_polym(newdata, degree = max(attr(object, "degree")),
		   vals_per_dim = attr(object, "vals_per_dim"),
	     coefs = attr(object, "coefs"))
    }
}
