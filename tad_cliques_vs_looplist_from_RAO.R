overlap_corner_peaks <- read.table("IMR90_tads_intersect_with_loops_from_RAO.bed")

colnames(overlap_corner_peaks) <- c("chr","start","end","clique_size","tad_id",
                                    "overlap")

singleton <- overlap_corner_peaks[overlap_corner_peaks$clique_size==1,]
binary <- overlap_corner_peaks[overlap_corner_peaks$clique_size==2,]
geq2 <- overlap_corner_peaks[overlap_corner_peaks$clique_size>=2,]
geq3 <- overlap_corner_peaks[overlap_corner_peaks$clique_size>=3,]
geq4 <- overlap_corner_peaks[overlap_corner_peaks$clique_size>=4,]
geq5 <- overlap_corner_peaks[overlap_corner_peaks$clique_size>=5,]
geq6 <- overlap_corner_peaks[overlap_corner_peaks$clique_size>=6,]
geq7 <- overlap_corner_peaks[overlap_corner_peaks$clique_size>=7,]
geq8 <- overlap_corner_peaks[overlap_corner_peaks$clique_size>=8,]
geq9 <- overlap_corner_peaks[overlap_corner_peaks$clique_size>=9,]
geq10 <- overlap_corner_peaks[overlap_corner_peaks$clique_size>=10,]
geq11 <- overlap_corner_peaks[overlap_corner_peaks$clique_size>=11,]
geq12 <- overlap_corner_peaks[overlap_corner_peaks$clique_size>=12,]
geq13 <- overlap_corner_peaks[overlap_corner_peaks$clique_size>=13,]
  
overlap_cat <- data.frame(clique_size=c("1","==2",">=2",">=3",">=4",
                                        ">=5",">=6",">=7",">=8",">=9",
                                        ">=10",">=11",">=12",">=13"),
                          overlapping_count=c(1:14),total_count=c(1:14))

overlap_cat[1,"overlapping_count"] <- table(singleton$overlap)[2]
overlap_cat[1,"total_count"] <- nrow(singleton)
overlap_cat[2,"overlapping_count"] <- table(binary$overlap)[2]
overlap_cat[2,"total_count"] <- nrow(binary)

overlap_cat[overlap_cat$clique_size==">=2","overlapping_count"] <- table(geq2$overlap)[2]
overlap_cat[overlap_cat$clique_size==">=2","total_count"] <- nrow(geq2)

overlap_cat[overlap_cat$clique_size==">=3","overlapping_count"] <- table(geq3$overlap)[2]
overlap_cat[overlap_cat$clique_size==">=3","total_count"] <- nrow(geq3)

overlap_cat[overlap_cat$clique_size==">=4","overlapping_count"] <- table(geq4$overlap)[2]
overlap_cat[overlap_cat$clique_size==">=4","total_count"] <- nrow(geq4)

overlap_cat[overlap_cat$clique_size==">=5","overlapping_count"] <- table(geq5$overlap)[2]
overlap_cat[overlap_cat$clique_size==">=5","total_count"] <- nrow(geq5)

overlap_cat[overlap_cat$clique_size==">=6","overlapping_count"] <- table(geq6$overlap)[2]
overlap_cat[overlap_cat$clique_size==">=6","total_count"] <- nrow(geq6)

overlap_cat[overlap_cat$clique_size==">=7","overlapping_count"] <- table(geq7$overlap)[2]
overlap_cat[overlap_cat$clique_size==">=7","total_count"] <- nrow(geq7)

overlap_cat[overlap_cat$clique_size==">=8","overlapping_count"] <- table(geq8$overlap)[2]
overlap_cat[overlap_cat$clique_size==">=8","total_count"] <- nrow(geq8)

overlap_cat[overlap_cat$clique_size==">=9","overlapping_count"] <- table(geq9$overlap)[2]
overlap_cat[overlap_cat$clique_size==">=9","total_count"] <- nrow(geq9)

overlap_cat[overlap_cat$clique_size==">=10","overlapping_count"] <- table(geq10$overlap)[2]
overlap_cat[overlap_cat$clique_size==">=10","total_count"] <- nrow(geq10)

overlap_cat[overlap_cat$clique_size==">=11","overlapping_count"] <- table(geq11$overlap)[2]
overlap_cat[overlap_cat$clique_size==">=11","total_count"] <- nrow(geq11)

overlap_cat[overlap_cat$clique_size==">=12","overlapping_count"] <- table(geq12$overlap)[2]
overlap_cat[overlap_cat$clique_size==">=12","total_count"] <- nrow(geq12)

overlap_cat[overlap_cat$clique_size==">=13","overlapping_count"] <- table(geq13$overlap)[2]
overlap_cat[overlap_cat$clique_size==">=13","total_count"] <- nrow(geq13)

overlap_cat$percent <- (overlap_cat$overlapping_count/overlap_cat$total_count)*100

overlap_cat$clique_size <- factor(overlap_cat$clique_size, 
                                 levels = c("1","==2",">=2",">=3",">=4",
                                            ">=5",">=6",">=7",">=8",">=9",
                                            ">=10",">=11",">=12",">=13"))

ggplot(data = overlap_cat, aes(x=clique_size,y=percent)) + 
  geom_point(color="darkgreen",size=3) + 
  ylim(10,30) +
  geom_vline(xintercept = 3.5) +
  theme_classic() + 
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 20)) 
