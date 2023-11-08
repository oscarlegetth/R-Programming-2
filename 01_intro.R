pbmc.data = readRDS("data/pbmc.rds")
dim(pbmc.data)

class(pbmc.data)
pbmc.data = as.matrix(pbmc.data)

# sum along each column to get the number of molecules found in each cell
n_count_rna = apply(pbmc.data, MARGIN = 2, sum)
# n_count_rna

# helper function to count the number of elements in an array that are above zero
count_above_zero = function(arr)
{
  n_above_zero = 0
  for (x in arr)
  {
    if (x > 0)
    {
      n_above_zero = n_above_zero + 1
    }
  }
  return(n_above_zero)
}

#n_feature_rna = apply(pbmc.data, MARGIN = 2, count_above_zero)
# better solution:
n_feature_rna = apply(pbmc.data > 0, MARGIN = 2, sum)

# calculate number of counts that where from mitochondrial genes
# grep("^MT-", row.names((pbmc.data))) returns row names that start with "MT-"
n_mt_count_rna = apply(pbmc.data[grep("^MT-", row.names((pbmc.data))),], MARGIN = 2, sum)
# element wise divide total mitochondrial counts with total counts to get our percentage
mt_percent_rna = n_mt_count_rna / n_count_rna


library(ggplot2)
library(gridExtra)
df = data.frame(n_features = n_feature_rna, n_count = n_count_rna, mt_percent = mt_percent_rna)

# make violin plots
plots = lapply(names(df), function(category)
{
  ggplot(df, aes(x=1, y=df[[category]])) + 
  geom_violin() +
  geom_jitter(shape=".", position=position_jitter(0.2)) +
  labs(x = category, y = "the rest of the plot")
})

grid.arrange(grobs = plots, ncol = length(names(df)))

# make scatter plots of n_features and mt_percent vs n_count
plots = lapply(c("n_features", "mt_percent"), function(category)
{
  ggplot(df, aes(x=df[["n_count"]], y=df[[category]])) + 
  geom_point() + 
  labs(x = "n_count", y = category)
})

grid.arrange(grobs = plots, ncol = 2)

# filter bad rows and columns
library(dplyr)
n_cells_per_feature = apply(pbmc.data > 0, 1, sum)
df.filtered = df %>% filter(n_features >= 200 & n_features <= 2500 & mt_percent < 0.05)
pbmc.filtered = pbmc.data[names(which(n_cells_per_feature > 2)), rownames(df.filtered)]
dim(pbmc.filtered)
