# Main
library(tidyverse)
library(Matrix)
library(readxl)
library(FactoMineR)
library(factoextra)
library(rstatix)

# Seq data
# library(DropletUtils) # Bioconductor package
library(Seurat)
library(DESeq2) # Bioconductor package

# Plot types
library(kableExtra)
library(platetools)
library(ComplexHeatmap) # Bioconductor package
library(vsn) # Bioconductor package
library(dendsort)
library(ggplotify)
library(umap)
library(corrplot)
library(ggcorrplot)
library(plotly)
library(rgl)
# library(VennDiagram)
library(ggvenn)
library(gt)
library(mmtable2)

# Colours
library(viridis)
library(scico)
library(RColorBrewer)

# Plot aestetics
library(ggthemes)
library(scales)
library(ggrepel)
library(ggpubr)
library(greekLetters)
library(patchwork)
library(cowplot)
library(ggh4x)
library(geomtextpath)
library(ggbreak)
library(MASS)

# Other
library(apeglm) # Bioconductor package
library(ashr)


## From online; Slightly modified from BUSpaRse, just to avoid installing a few dependencies not used here
read_count_output <- function(dir, name) {
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/", name, ".mtx"))
  m <- Matrix::t(m)
  m <- as(m, "dgCMatrix")
  # The matrix read has cells in rows
  ge <- ".genes.txt"
  genes <- readLines(file(paste0(dir, "/", name, ge)))
  barcodes <- readLines(file(paste0(dir, "/", name, ".barcodes.txt")))
  colnames(m) <- barcodes
  rownames(m) <- genes
  return(m)
}

## For re-ordered antibody plot
scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}
reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

# From online, for using a table in a combined figure
gt_grob <- function(gt_object, ...){

  out_name <- file.path(
    tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")
  )

  gtsave(gt_object, out_name, ...)

  # in_png <- png::readPNG(out_name)
  in_png <- magick::image_read(out_name) %>%
    magick::image_ggplot(interpolate = TRUE)

  on.exit(file.remove(out_name), add=TRUE)

  # grid::rasterGrob(in_png)
  return(in_png)

}
