#' coxpf: Profile Monitoring Methods for Cox Modeling
#'
#' Functions to implement nonparametric profile monitoring methods for the
#' cumulative baseline hazard for Cox models. 
#' Major quantities that are computed:
#' \itemize{
#' \item reference shape when a matrix of shapes are given
#' \item reference variabilities after the reference shape is subtracted
#' from all shapes
#' }
#'
#' @docType package
#' @name coxpf
#' @import survival
#' @useDynLib coxpf
NULL