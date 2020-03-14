clic_file <- read.table("IMR90_clicum.tadid")
orient <- read.table("TADs_with_boundary_orientation_no_cen.csv")
clic_orient <- merge(clic_file,orient,by.x = "V5",by.y = "V1")
colnames(clic_orient) <- c("tadid","chr","start","end","cliq_size","orient")


convergents <- clic_orient[clic_orient$orient=="F,R",c("chr","start","end","cliq_size","orient")]
colnames(convergents) <- NULL
write.table(convergents,"IMR90_convergents.tsv",quote = F,sep = "\t",row.names = F,col.names = F)

#Wrote AWK after this 
