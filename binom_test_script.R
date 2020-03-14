binom_test_for_cliq_against_global <- function(clic_orient)
{
  stats <- data.frame(cliq_cat=c(1:10),confident1=c(1:10),confident2=c(1:10),convergent_tads=c(1:10),
                      total_tads_in_cat=c(1:10),global_convergent_tads=c(1:10),pvalue=c(1:10))
  rownames(stats) <- c("==2",">=2",">=3",">=4",">=5",">=6",">=7",">=8",">=9",">=10")
  global_avg <- nrow(clic_orient[clic_orient$orient=="F,R",])/nrow(clic_orient)
  
  k <- binom.test(nrow(clic_orient[clic_orient$cliq_size==2 & clic_orient$orient=="F,R",]),
                  nrow(clic_orient[clic_orient$cliq_size==2,]),p=global_avg)
  stats["==2","cliq_cat"] <- "==2"
  stats["==2","confident1"] <- k$conf.int[1]
  stats["==2","confident2"] <- k$conf.int[2]
  stats["==2","convergent_tads"] <- k$statistic
  stats["==2","total_tads_in_cat"] <- k$parameter
  stats["==2","global_convergent_tads"] <- k$null.value
  stats["==2","pvalue"] <- k$p.value
  
  for (i in c(2:10)) 
  {
    paste(">=",i,sep = "")
    k <- binom.test(nrow(clic_orient[clic_orient$cliq_size>=i & clic_orient$orient=="F,R",]),
                    nrow(clic_orient[clic_orient$cliq_size>=i,]),p=global_avg)
    stats[paste(">=",i,sep = ""),"cliq_cat"] <- paste(">=",i,sep = "")
    stats[paste(">=",i,sep = ""),"confident1"] <- k$conf.int[1]
    stats[paste(">=",i,sep = ""),"confident2"] <- k$conf.int[2]
    stats[paste(">=",i,sep = ""),"convergent_tads"] <- k$statistic
    stats[paste(">=",i,sep = ""),"total_tads_in_cat"] <- k$parameter
    stats[paste(">=",i,sep = ""),"global_convergent_tads"] <- k$null.value
    stats[paste(">=",i,sep = ""),"pvalue"] <- k$p.value
  }
  return(stats)
}

run_binom_file <- function(tadid_file,orient_file,outputfile)
{
tadfile <- read.table(tadid_file)
orientfile <- read.table(orient_file)
clic_orient <- merge(tadfile,orientfile,by.x = "V5",by.y = "V1")
colnames(clic_orient) <- c("tadid","chr","start","end","cliq_size","orient")
stats_df <- binom_test_for_cliq_against_global(clic_orient)
write.table(x = stats_df,file = outputfile,quote = F,sep = "\t",row.names = F, col.names = T)
}

run_binom_file(tadid_file = "IMR90_clicum.tadid",orient_file = "TADs_with_boundary_orientation_no_cen.csv",outputfile = "IMR90_binom.tsv")
run_binom_file(tadid_file = "HMEC_clicnum.tadid",orient_file = "HMEC_TADs_with_boundary_orientation_no_cen.csv",outputfile = "HMEC_binom.tsv")
run_binom_file(tadid_file = "HUVEC_clicnum.tadid",orient_file = "HUVEC_TADs_with_boundary_orientation_no_cen.csv",outputfile = "HUVEC_binom.tsv")
run_binom_file(tadid_file = "K562_clicnum.tadid",orient_file = "K562_TADs_with_boundary_orientation_no_cen.csv",outputfile = "K562_binom.tsv")


