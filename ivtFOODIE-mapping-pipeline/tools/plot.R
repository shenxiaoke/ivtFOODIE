library(ggplot2)
library(reshape2)

# Optional CLI argument: path to info_merge directory.
# Default (no argument) works with Snakefile rule that cds into info_merge.
args <- commandArgs(trailingOnly = TRUE)
info_dir <- if (length(args) >= 1 && nchar(args[1]) > 0) args[1] else "."
info_dir <- normalizePath(info_dir)

datasets <- list.dirs(path = info_dir, recursive = FALSE, full.names = TRUE)
# list.dirs includes the root path; drop it so we only iterate samples.
datasets <- datasets[datasets != info_dir]

pdf(file.path(info_dir, "converted_ratio.pdf"), width = 3.7, height = 2.5)
print(datasets)

for (f in datasets) {
  prefix <- basename(f)
  fname <- file.path(f, paste0(prefix, ".call.o.txt"))
  print(fname)

  if (!file.exists(fname)) {
    next
  }

  info <- file.info(fname)
  empty <- (info$size[1] == 0)

  if (empty) {
    next
  }

  print(fname)

  stats <- read.delim(fname, header = FALSE)
  # stats <- as.numeric(unlist(stats[2,]))
  stats <- as.numeric(unlist(stats[1, ]))
  stats <- matrix(stats, nrow = 8)
  stats <- t(stats[c(2, 4, 6, 8), ] / (stats[c(1, 3, 5, 7), ] + stats[c(2, 4, 6, 8), ]))
  print(stats)
  rownames(stats) <- colnames(stats) <- c("A", "T", "C", "G")
  stats <- as.data.frame(stats)
  stats$preceding_base <- rownames(stats)

  # prefix <- substring(prefix,8)

  stats <- melt(stats, variable.name = "following_base", value.name = "converted_ratio")
  p1 <- ggplot(data = stats, aes(x = following_base, y = preceding_base)) +
    geom_tile(aes(fill = converted_ratio)) +
    geom_text(aes(label = round(converted_ratio, 3))) +
    theme_classic() +
    theme(axis.line = element_blank(), plot.title = element_text(size = 5, face = "bold")) +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) +
    ggtitle(prefix)
  print(p1)
}

dev.off()
