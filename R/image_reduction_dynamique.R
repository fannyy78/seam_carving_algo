#' Calcul de l'energie dynamique d'une image
#'
#' @param gray_img Image en niveaux de gris (cimg)
#' @return Image d'energie (cimg)
#' @export
compute_dynamic_energy <- function(gray_img) {
  gx <- imgradient(gray_img, "x")
  gy <- imgradient(gray_img, "y")
  energy <- sqrt(gx^2 + gy^2)
  return(energy)
}

#' Calcule la matrice cumulative d'energie dynamique
#'
#' @param energy_mat Matrice d'energie (h x w)
#' @return Matrice cumulative (h x w)
#' @export
cumulative_dynamic_energy <- function(energy_mat) {
  h <- nrow(energy_mat)
  w <- ncol(energy_mat)
  cum_energy <- energy_mat

  for (i in 2:h) {
    for (j in 1:w) {
      left   <- if (j > 1) cum_energy[i - 1, j - 1] else Inf
      middle <- cum_energy[i - 1, j]
      right  <- if (j < w) cum_energy[i - 1, j + 1] else Inf
      cum_energy[i, j] <- energy_mat[i, j] + min(left, middle, right)
    }
  }
  return(cum_energy)
}

#' Retrouve la seam minimale dynamique
#'
#' @param cum_energy Matrice cumulative d'energie (h x w)
#' @return Vecteur (longueur h) representant la seam
#' @export
find_dynamic_seam <- function(cum_energy) {
  h <- nrow(cum_energy)
  w <- ncol(cum_energy)
  seam <- integer(h)
  seam[h] <- which.min(cum_energy[h, ])

  for (i in (h - 1):1) {
    prev_col <- seam[i + 1]
    options <- c(
      if (prev_col > 1) cum_energy[i, prev_col - 1] else Inf,
      cum_energy[i, prev_col],
      if (prev_col < w) cum_energy[i, prev_col + 1] else Inf
    )
    offset <- which.min(options) - 2
    seam[i] <- prev_col + offset
  }
  return(seam)
}

#' Supprime une couture verticale dynamique
#'
#' @param im Image (objet `cimg`)
#' @param seam Vecteur des indices de colonnes a supprimer
#' @return Image modifiee (objet `cimg`)
#' @export
remove_dynamic_vertical_seam <- function(im, seam) {
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

#' Reduction d'image par seam carving dynamique
#'
#' @param im Image a reduire (objet `cimg`)
#' @param n_seams Nombre de coutures a supprimer
#'
#' @return Image reduite (objet `cimg`)
#' @export
reduce_dynamic_image_seams <- function(im, n_seams) {
  par(mfrow = c(1, 2))  # Affichage avant/apres
  plot(im, main = "Image originale")

  for (i in 1:n_seams) {
    cat("Suppression seam dynamique", i, "\n")
    gray <- grayscale(im)
    energy <- compute_dynamic_energy(gray)
    energy_mat <- t(as.matrix(energy))

    cum_energy <- cumulative_dynamic_energy(energy_mat)
    best_seam <- find_dynamic_seam(cum_energy)
    im <- remove_dynamic_vertical_seam(im, best_seam)
  }

  plot(im, main = paste("Image apres reduction -", n_seams, "seams supprimees"))

  return(im)
}
