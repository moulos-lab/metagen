library(rmarkdown)
library(ggplot2)
library(dada2)
library(phyloseq)
library(Biostrings)
library(kableExtra)

ps <- readRDS("ps.rds")
enzymes <- read.delim("~/Documents/run343_16s/picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv")
pathways <- read.delim("~/Documents/run343_16s/picrust2_out_pipeline/pathways_out/path_abun_unstrat_descrip.tsv")

#1
# Alpha Diversity
plot_richness(
	ps, 
	color="Subject",
	measures=c("Chao1", "Shannon"),
	title = "Alpha diversity"
	) + 
  geom_point(size=5, alpha=0.7)

#2
# Plot Beta diversity (Ordination)
# OTUs apearing at least 3 times in at least half of the samples
# Plot only top-10 phyla of these OTUs
wh0 = genefilter_sample(
	ps, 
	filterfun_sample(function(x) x > 3), 
	A=0.5*nsamples(ps)
	)
GP1 = prune_taxa(wh0, ps)
GP1 = ps
GP1 = transform_sample_counts(GP1, function(x) 1E6 * x/sum(x))
phylum.sum = tapply(
	taxa_sums(GP1), 
	tax_table(GP1)[, "Order"], 
	sum, 
	na.rm=TRUE
	)
top10phyla = names(sort(phylum.sum, TRUE))[1:10]
GP1 = prune_taxa((tax_table(GP1)[, "Order"] %in% top10phyla), GP1)
GP.ord <- ordinate(GP1, "MDS", "bray")

plot_ordination(
	GP1, 
	GP.ord, 
	type="split", 
	color="Order", 
	label="Subject",
	title="Beta diversity"
	) + geom_point(size=2) + scale_color_brewer(palette="Spectral")

#3
# Abundance barplot based in top30 appearing OTUs accross all samples
top30 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:30]
ps_top30 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps_top30 <- prune_taxa(top30, ps_top30)

plot_bar(
	ps_top30, 
	"Family", 
	fill="Genus", 
	facet_grid=~Subject,
	title="Abundance barplot"
	) +
scale_fill_brewer(palette="Paired")

#4
# Tree
# Top 80 listed OTUs
ps_tree = prune_taxa(taxa_names(ps)[1:80], ps)

plot_tree(
  ps_tree,
  ladderize="left",
  shape="Sample",
  color="Phylum",
  label.tips = "Class",
  size="abundance",
  base.spacing = 0.04,
  method = "sampledodge",
  sizebase = 10,
  title = "Phyla-Sample Clustering"
  ) + coord_polar(theta="y")

#5
# Net
# top-100 appearing OTUs by phylum
ps_net <- prune_taxa(names(sort(taxa_sums(ps), TRUE))[1:100], ps)
jg <- make_network(GP100, "taxa", "jaccard", 0.3)
plot_network(title = "Bacterial classes co-occurence network", jg, GP100, "taxa", color = "Phylum", line_weight = 0.4, label = "Class")

#6
# Enzymes Treemap
class <- case_when(
  grepl(pattern = "EC:1", x = enzymes$function.) ~ "Oxidoreductases",
  grepl(pattern = "EC:2", x = enzymes$function.) ~ "Transferases",
  grepl(pattern = "EC:3", x = enzymes$function.) ~ "Hydrolases",
  grepl(pattern = "EC:4", x = enzymes$function.) ~ "Lyases",
  grepl(pattern = "EC:5", x = enzymes$function.) ~ "Isomerases",
  grepl(pattern = "EC:6", x = enzymes$function.) ~ "Ligases",
  grepl(pattern = "EC:7", x = enzymes$function.) ~ "Translocases"
)

enzymes <- add_column(enzymes, class, .after = 2)
enzymes_expanded <- gather(enzymes, "sample", "abundance", 4:ncol(enzymes))

treemap(enzymes_expanded,
        aspRatio = 2,
        index = c("sample",
                  "class"),
        vSize ="abundance", 
        vColor = "sample",
        type="categorical",
        align.labels = list(c("centre",
                              "bottom"),
                            c("left",
                              "top")),
        algorithm = "pivotSize" ,
        sortID = "abundance", 
        mirror.x = T,
        palette = "Paired",
        title = "Enzyme class abundance",
        fontcolor.labels = c("#435856","white"), 
        fontsize.labels = c(15,9),
        fontface.labels = 2,
        border.col = "grey", 
        bg.labels = 220)

#7
# Pathway Heatmap
pt <- pathways
pt[] <- paste("description:",pathways$description)
rownames(pathways) <- pathways[,1]
heatmaply(normalize(pathways[,3:ncol(pathways)]), 
          custom_hovertext = pt, 
          xlab = "Samples", 
          ylab = "KEGG Pathway", 
          main = "Pathways abundance")

#8
# Data trunc
kable(data_trunc) %>%
kable_styling("striped", full_width = F) %>%
column_spec(5, color = "green", bold = TRUE)