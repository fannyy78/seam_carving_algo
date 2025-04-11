#' Calcul de l'energie d'une image
#' @param gray_img Image en niveaux de gris (cimg)
#' @return Image d'energie
#'@export
compute_energy <- function(gray_img) {
  gx <- imgradient(gray_img, "x")
  gy <- imgradient(gray_img, "y")
  energy <- sqrt(gx^2 + gy^2)
  return(energy)
}

#' Trace une couture gloutonne
#' @param energy_mat Matrice d'energie (h x w)
#' @param start_col Colonne de depart
#' @param penalty Penalite pour diagonales
#' @return Vecteur d'indices de colonnes
#' @export
trace_greedy_seam <- function(energy_mat, start_col, penalty = 0) {
  h <- nrow(energy_mat)
  w <- ncol(energy_mat)
  seam <- integer(h)
  seam[1] <- start_col

  for (i in 2:h) {
    prev_col <- seam[i - 1]
    values <- c(
      if (prev_col > 1) energy_mat[i, prev_col - 1] + penalty else Inf,
      energy_mat[i, prev_col],
      if (prev_col < w) energy_mat[i, prev_col + 1] + penalty else Inf
    )
    offset <- which.min(values) - 2
    seam[i] <- prev_col + offset
  }
  return(seam)
}

#' Supprime une couture verticale
#' @param im Image (objet `cimg`)
#' @param seam Vecteur d'indices de colonne a supprimer
#'
#' @return Image (objet `cimg`) avec 1 colonne en moins
#' @export
remove_vertical_seam <- function(im, seam) {
  im_array <- as.array(im)
  w <- dim(im_array)[1]
  h <- dim(im_array)[2]
  c <- dim(im_array)[4]

  new_array <- array(0, dim = c(w - 1, h, 1, c))

  for (y in 1:h) {
    col <- seam[y]
    keep_cols <- setdiff(1:w, col)
    for (ch in 1:c) {
      new_array[, y, 1, ch] <- im_array[keep_cols, y, 1, ch]
    }
  }

  return(as.cimg(new_array))
}

#' Applique le seam carving
#' @param im Image (objet `cimg`)
#' @param n_seams Nombre de coutures a supprimer
#'
#' @return Image reduite
#' @export
seam_carving <- function(im, n_seams = 10) {
  for (step in 1:n_seams) {
    cat("Suppression seam", step, "\n")
    gray <- grayscale(im)
    energy <- compute_energy(gray)
    energy_mat <- t(as.matrix(energy))

    h <- nrow(energy_mat)
    w <- ncol(energy_mat)

    all_seams <- vector("list", w)
    for (col in 1:w) {
      all_seams[[col]] <- trace_greedy_seam(energy_mat, col)
    }

    seam_energies <- sapply(all_seams, function(seam) {
      sum(mapply(function(i, j) energy_mat[i, j], 1:h, seam))
    })
    best_idx <- which.min(seam_energies)
    best_seam <- all_seams[[best_idx]]

    im <- remove_vertical_seam(im, best_seam)
  }
  return(im)
}

#' Reduit une image avec seam carving
#'
#' @param im Image (objet `cimg`)
#' @param n_seams Nombre de coutures a supprimer
#'
#' @return Image reduite
#' @export
reduce_image_seams <- function(im, n_seams = 10) {
  cat("Dimensions AVANT reduction :", dim(im), "\n")
  im_reduite <- seam_carving(im, n_seams)
  cat("Dimensions APRES reduction :", dim(im_reduite), "\n")

  par(mfrow = c(1, 2))
  plot(im, main = "Image originale")
  plot(im_reduite, main = paste("Image reduite (", n_seams, "seams)", sep = ""))

  return(im_reduite)
}

#' @useDynLib M2algorithmique
#' @importFrom Rcpp sourceCpp
#' @importFrom graphics par
NULL
