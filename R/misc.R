# #' Check whether the input is a valid Hi-C
# #'
# check_valid_hic <- function(hic_matrix) {
#   chrom <- attr(hic_matrix, "chrom")
#   resol <- attr(hic_matrix, "resol")
#   type <- attr(hic_matrix, "type")
# 
#   if (length(chrom) < 1)
#     stop("chrom is missing")
# 
#   if (is.null(resol))
#     stop("resol is missing")
#   if (length(resol) != 1 ||
#       !(resol %in% supported_resol)) {
#     stop(str_interp("Invalid resolution: ${resol}"))
#   }
# 
#   if (length(type) != 1)
#     stop("Invalid matrix type")
# 
#   stopifnot(colnames(hic_matrix)[1:5] == c("chrom1", "pos1", "chrom2", "pos2", "score"))
# 
#   TRUE
# }
