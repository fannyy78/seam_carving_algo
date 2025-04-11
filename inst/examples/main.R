devtools::install_github("vrunge/M2algorithmique")

# Charge les libs
library(Rcpp)
library(imager)
library(M2algorithmique)

ls("package:M2algorithmique")

setwd("C:/Alix/MASTER/Algorithmique/M2algorithmique")
devtools::document()
devtools::load_all()
tools::showNonASCIIfile("R/image_reduction.R")


devtools::check()

im <- load.example("parrots")

#R  naif
im_reduite <- reduce_image_seams(im, 5)

#C++ naif
im_reduite <- greedy_carving_cpp(im, 30)
plot(im_reduite)

#R_Dynamique
n_seams <- 100
im_reduite <- reduce_dynamic_image_seams(im, n_seams)

#C++_Dynamique
im_reduite <- seam_carving_dp(im, 100)
plot(im_reduite, main = "Image réduite (Programmation Dynamique C++)")

















