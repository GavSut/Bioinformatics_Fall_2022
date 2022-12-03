#Ordinations

#Ordinate
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

#Plotting Ordination
plot_ordination(ps.prop, ord.nmds.bray, color="Spider", title="Bray NMDS")+
  theme_classic()

ord.pcoa.bray <- ordinate(ps.prop, method="PCoA", distance="bray")
ord.pcoa.bray
plot_ordination(ps.prop, ord.pcoa.bray, type = "taxa", color="Spider", title="Bray PCoa")+
  theme_classic()


GP = ps
wh0 = genefilter_sample(GP, filterfun_sample(function(x) x > 5), A=0.5*nsamples(GP))
GP1 = prune_taxa(wh0, GP)

GP1 = transform_sample_counts(GP1, function(x) 1E6 * x/sum(x))

phylum.sum = tapply(taxa_sums(GP1), tax_table(GP1)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
GP1 = prune_taxa((tax_table(GP1)[, "Phylum"] %in% top5phyla), GP1)

GP.ord <- ordinate(GP1, "PCoa", "bray")
p1 = plot_ordination(GP1, GP.ord, type="taxa", color="Phylum", title="taxa")
p1  
p1 + facet_wrap(~Phylum, 3)

p2 = plot_ordination(GP1, GP.ord, type="samples", color="Spider") 
p2

p3 = plot_ordination(GP1, GP.ord, type="biplot", color="samples", title="biplot")
GP1.shape.names = get_taxa_unique(GP1, "Phylum")
GP1.shape <- 15:(15 + length(GP1.shape.names) - 1)
names(GP1.shape) <- GP1.shape.names
GP1.shape["samples"] <- 16
p3 + scale_shape_manual(values=GP1.shape)











