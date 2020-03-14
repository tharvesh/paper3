clic_file <- read.table("IMR90_clicum.tadid")
orient <- read.table("TADs_with_boundary_orientation_no_cen.csv")
clic_orient <- merge(clic_file,orient,by.x = "V5",by.y = "V1")
colnames(clic_orient) <- c("tadid","chr","start","end","cliq_size","orient")

binary_interact_all_ctcf <- clic_orient[clic_orient$cliq_size==2,] 

binary_interact_bed <- read.table("binary_interactions_dist.bed")
colnames(binary_interact_bed) <- c("chrA","s1","e1","chrB","s2","e2","dist","chrC",
                                   "s3","e3","clique_size","tadid")


binary_interact_bed_with_orient <- merge(x = binary_interact_bed,y = clic_orient,by = "tadid")

bin_0_500kb <- nrow(binary_interact_bed_with_orient[binary_interact_bed_with_orient$dist<500000 & binary_interact_bed_with_orient$orient=="F,R",])/
  nrow(binary_interact_bed_with_orient[binary_interact_bed_with_orient$dist<500000,])*100
bin_500kb_1mb <- nrow(binary_interact_bed_with_orient[binary_interact_bed_with_orient$dist>=500000 & binary_interact_bed_with_orient$dist<1000000 & binary_interact_bed_with_orient$orient=="F,R",])/
  nrow(binary_interact_bed_with_orient[binary_interact_bed_with_orient$dist>=500000 & binary_interact_bed_with_orient$dist<1000000,])*100
bin_1mb_2mb <- nrow(binary_interact_bed_with_orient[binary_interact_bed_with_orient$dist>=1000000 & binary_interact_bed_with_orient$dist<2000000 & binary_interact_bed_with_orient$orient=="F,R",])/
  nrow(binary_interact_bed_with_orient[binary_interact_bed_with_orient$dist>=1000000 & binary_interact_bed_with_orient$dist<2000000,])*100
bin_gr_2mb <- nrow(binary_interact_bed_with_orient[binary_interact_bed_with_orient$dist>2000000 & binary_interact_bed_with_orient$orient=="F,R",])/
  nrow(binary_interact_bed_with_orient[binary_interact_bed_with_orient$dist>2000000,])*100

bins <- c("0-500Kb","500Kb-1Mb","1Mb-2Mb",">=2Mb")
bin_val <- c(bin_0_500kb,bin_500kb_1mb,bin_1mb_2mb,bin_gr_2mb)

plot_df <- data.frame(bins=bins,bin_val=bin_val)

plot_df$bins <- factor(plot_df$bins, levels = c("0-500Kb","500Kb-1Mb","1Mb-2Mb",">=2Mb"))

ggplot(data = plot_df, aes(x=bins,y=bin_val,group=1)) + 
  xlab("Genomic distance between binary interactions") +
  ylab("Percent of convergent motifs") +
  geom_line(linetype = "dashed",color="darkblue",size=1) +
  geom_point(color="darkblue",size=3) + 
  ylim(5,30) +
  theme_classic() + 
  theme(axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 14,angle = 90,hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 15)) 

bin_gr_500kb <- nrow(binary_interact_bed_with_orient[binary_interact_bed_with_orient$dist>=500000 & binary_interact_bed_with_orient$orient=="F,R",])/
  nrow(binary_interact_bed_with_orient[binary_interact_bed_with_orient$dist>=500000,])*100
bar_df <- data.frame(bins=c(bins[1],">=500Kb"),bin_val=c(bin_val[1],bin_gr_500kb))

bar_df$bins <- factor(bar_df$bins,levels = c(bins[1],">=500Kb"))

ggplot(data = bar_df, aes(x=bins,y = bin_val)) + geom_bar(stat="identity", fill="steelblue") +
  xlab("Genomic distance between binary interactions") +
  ylab("Percent of convergent motifs") +
  theme_minimal() + 
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20))

ggsave("percent_conv_motif_in_binary_interaction.pdf", plot = last_plot(), device = "pdf",
       scale = 1, width = 7, height = 7, units = "in",
       dpi = 300, limitsize = TRUE)

bin_int_conv_motif_le_500kb <- nrow(binary_interact_bed_with_orient[binary_interact_bed_with_orient$dist<500000 & binary_interact_bed_with_orient$orient=="F,R",])
not_bin_int_conv_motif_le_500kb <- nrow(binary_interact_bed_with_orient[binary_interact_bed_with_orient$dist<500000,])-bin_int_conv_motif_le_500kb
bin_int_conv_motif_geq_500kb <- nrow(binary_interact_bed_with_orient[binary_interact_bed_with_orient$dist>=500000 & binary_interact_bed_with_orient$orient=="F,R",])
not_bin_int_conv_motif_geq_500kb <- nrow(binary_interact_bed_with_orient[binary_interact_bed_with_orient$dist>=500000,])-bin_int_conv_motif_geq_500kb

fischer_mat <- matrix(c(bin_int_conv_motif_le_500kb,bin_int_conv_motif_geq_500kb,not_bin_int_conv_motif_le_500kb,not_bin_int_conv_motif_geq_500kb),nrow = 2, byrow = 1)
fisher.test(fischer_mat)




