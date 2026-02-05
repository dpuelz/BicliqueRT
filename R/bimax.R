# Internal Bimax biclustering functions
# Replaces the biclust package dependency with an embedded implementation
# of the Bimax algorithm (Prelic et al., 2006)

#' Run Bimax biclustering on a binary matrix
#'
#' @param x A binary matrix (0/1 entries).
#' @param minr Minimum number of rows in a bicluster.
#' @param minc Minimum number of columns in a bicluster.
#' @param number Maximum number of biclusters to find.
#'
#' @return A list with components:
#'   \item{RowxNumber}{Logical matrix (nrow x Number), rows belonging to each bicluster}
#'   \item{NumberxCol}{Logical matrix (Number x ncol), cols belonging to each bicluster}
#'   \item{Number}{Integer, number of biclusters found}
#'
#' @keywords internal
bimax_biclust <- function(x, minr, minc, number = 100) {
  x <- as.matrix(x)
  nr <- nrow(x)
  nc <- ncol(x)

  # Ensure binary integer matrix
  storage.mode(x) <- "integer"

  result_rows <- integer(nr * number)
  result_cols <- integer(nc * number)
  num_found <- integer(1)

  out <- .C("bimax_c",
            data_vec = x,
            nrow_r = as.integer(nr),
            ncol_r = as.integer(nc),
            minr_r = as.integer(minr),
            minc_r = as.integer(minc),
            number_r = as.integer(number),
            result_rows = result_rows,
            result_cols = result_cols,
            num_found = num_found)

  nfound <- out$num_found

  if (nfound == 0) {
    RowxNumber <- matrix(FALSE, nrow = nr, ncol = 1)
    NumberxCol <- matrix(FALSE, nrow = 1, ncol = nc)
    return(list(RowxNumber = RowxNumber, NumberxCol = NumberxCol, Number = 0))
  }

  rmat <- matrix(out$result_rows, nrow = nr, ncol = number)[, 1:nfound, drop = FALSE]
  cmat <- matrix(out$result_cols, nrow = nc, ncol = number)[, 1:nfound, drop = FALSE]

  RowxNumber <- rmat == 1
  NumberxCol <- t(cmat == 1)

  rownames(RowxNumber) <- rownames(x)
  colnames(NumberxCol) <- colnames(x)

  list(RowxNumber = RowxNumber, NumberxCol = NumberxCol, Number = nfound)
}


#' Extract bicluster submatrices from a Bimax result
#'
#' Replaces \code{biclust::bicluster()}. Extracts the submatrix for each
#' bicluster found.
#'
#' @param x The original matrix that was biclustered.
#' @param bicResult A list returned by \code{bimax_biclust()}.
#'
#' @return A named list of submatrices, one per bicluster.
#'
#' @keywords internal
bimax_extract <- function(x, bicResult) {
  x <- as.matrix(x)
  nfound <- bicResult$Number
  result <- vector("list", nfound)

  for (k in seq_len(nfound)) {
    rows_k <- which(bicResult$RowxNumber[, k])
    cols_k <- which(bicResult$NumberxCol[k, ])
    if (length(rows_k) > 0 && length(cols_k) > 0) {
      result[[k]] <- x[rows_k, cols_k, drop = FALSE]
    } else {
      result[[k]] <- NULL
    }
  }

  names(result) <- paste0("Bicluster", seq_len(nfound))
  result
}


#' Extract row and column indices from a Bimax result
#'
#' Replaces \code{biclust::biclusternumber()}. Returns row and column indices
#' for each bicluster found.
#'
#' @param bicResult A list returned by \code{bimax_biclust()}.
#'
#' @return A named list where each element has \code{Rows} and \code{Cols}
#'   integer vectors.
#'
#' @keywords internal
bimax_indices <- function(bicResult) {
  nfound <- bicResult$Number
  result <- vector("list", nfound)

  for (k in seq_len(nfound)) {
    rows_k <- which(bicResult$RowxNumber[, k])
    cols_k <- which(bicResult$NumberxCol[k, ])
    result[[k]] <- list(Rows = rows_k, Cols = cols_k)
  }

  names(result) <- paste0("Bicluster", seq_len(nfound))
  result
}
