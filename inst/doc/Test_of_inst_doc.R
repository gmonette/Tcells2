#' ---
#' title: "Test to see what happens in inst/doc"
#' author: "Georges Monette"
#' date: "`r format(Sys.time(), '%H:%M %d %B %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      toc_depth: 4
#' ---
#' Generated:
{{format(Sys.time(), '%H:%M %d %B %Y')}}
#' 
#' <!--
#' Let's see whether this works to create a file in inst/doc
#' * update 2015 08 10: strarted this
#' -->
#' 
#+ opts,echo=FALSE
library(knitr)
opts_knit$set(width=120)
options(width=120)
opts_chunk$set(tidy=FALSE,comment='| ',fig.height=8,fig.width=10)
#'
data(cars)
head(cars)
plot( speed ~ dist, cars)