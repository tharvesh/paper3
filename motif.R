clic_file <- read.table("IMR90_clicum.tadid")
orient <- read.table("TADs_with_boundary_orientation_no_cen.csv")
clic_orient <- merge(clic_file,orient,by.x = "V5",by.y = "V1")
colnames(clic_orient) <- c("tadid","chr","start","end","cliq_size","orient")

cliq_ctcf_orient_table <- function(clic_orient_df,small_cliq_size,big_cliq_size,big_cliq_size_only_eq=FALSE) 
{
  small_clic_orient <- clic_orient_df[clic_orient_df$cliq_size==small_cliq_size,]
  if (big_cliq_size_only_eq==FALSE) {
    big_clic_orient <- clic_orient_df[clic_orient_df$cliq_size>=big_cliq_size,]
  }
  else
  {
    big_clic_orient <- clic_orient_df[clic_orient_df$cliq_size==big_cliq_size,]
  }
  
  small_clic_orient <- clic_orient_df[clic_orient_df$cliq_size==small_cliq_size,]
  big_clic_orient_table <- data.frame(table(big_clic_orient$orient))
  big_clic_orient_table$percent <- (big_clic_orient_table$Freq/sum(big_clic_orient_table$Freq))*100
  colnames(big_clic_orient_table) <- c("CategoryB","BigCliqFreq","BigCliqPercent")
  small_clic_orient_table <- data.frame(table(small_clic_orient$orient))
  small_clic_orient_table$percent <- (small_clic_orient_table$Freq/sum(small_clic_orient_table$Freq))*100
  colnames(small_clic_orient_table) <- c("Category","SmallCliqFreq","SmallCliqPercent")
  big_clic_orient_table$CategoryB <- NULL
  orient_percent_side <- cbind(small_clic_orient_table,big_clic_orient_table)
  
  return(orient_percent_side)
}

FR_orient_test <- function(orient_percent_side_table)
{
  ops <- orient_percent_side_table
  r1c1 <- ops[ops$Category=="F,R",2]
  r1c2 <- ops[ops$Category=="F,R",4]
  r2c1 <- sum(ops$SmallCliqFreq)-r1c1
  r2c2 <- sum(ops$BigCliqFreq)-r1c2
  
  mat <- matrix(c(r1c1,r1c2,r2c1,r2c2), nrow = 2, byrow = 1)
  test_result <- fisher.test(mat)
  
  res <- c(mat=c(r1c1,r1c2,r2c1,r2c2),smallCliq_FR_percent=ops[ops$Category=="F,R",3],BigCliq_FR_percent=ops[ops$Category=="F,R",5]
           ,statistics=test_result)
  
  return(res)
}

run_FR_orient_test <- function(clic_orient_df,small_cliq_size,big_cliq_size,big_cliq_size_only_eq=FALSE)
{
  ops <- cliq_ctcf_orient_table(clic_orient_df,small_cliq_size,big_cliq_size,big_cliq_size_only_eq)
  res <- FR_orient_test(ops)
}

fr_orient_result_1_gre2 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 1, big_cliq_size = 2)
fr_orient_result_1_eq2 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                             small_cliq_size = 1, big_cliq_size = 2,
                                             big_cliq_size_only_eq = TRUE)
fr_orient_result_1_gre3 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 1, big_cliq_size = 3)
fr_orient_result_1_gre4 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 1, big_cliq_size = 4)
fr_orient_result_1_gre5 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 1, big_cliq_size = 5)
fr_orient_result_1_gre6 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 1, big_cliq_size = 6)
fr_orient_result_1_gre7 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 1, big_cliq_size = 7)
fr_orient_result_1_gre8 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 1, big_cliq_size = 8)
fr_orient_result_1_gre9 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 1, big_cliq_size = 9)

fr_orient_result_1_gre10 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                               small_cliq_size = 1, big_cliq_size = 10)
fr_orient_result_1_gre11 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                               small_cliq_size = 1, big_cliq_size = 11)

fr_orient_result_1_gre12 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                               small_cliq_size = 1, big_cliq_size = 12)
fr_orient_result_1_gre13 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                               small_cliq_size = 1, big_cliq_size = 13)


fr_orient_result_2_gre2 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 2, big_cliq_size = 2)
fr_orient_result_2_gre3 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 2, big_cliq_size = 3)
fr_orient_result_2_gre4 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 2, big_cliq_size = 4)
fr_orient_result_2_gre5 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 2, big_cliq_size = 5)
fr_orient_result_2_gre6 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 2, big_cliq_size = 6)
fr_orient_result_2_gre7 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 2, big_cliq_size = 7)
fr_orient_result_2_gre8 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 2, big_cliq_size = 8)
fr_orient_result_2_gre9 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 2, big_cliq_size = 9)

fr_orient_result_2_gre10 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                               small_cliq_size = 2, big_cliq_size = 10)
fr_orient_result_2_gre11 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                               small_cliq_size = 2, big_cliq_size = 11)

fr_orient_result_2_gre12 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                               small_cliq_size = 2, big_cliq_size = 12)
fr_orient_result_2_gre13 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                               small_cliq_size = 2, big_cliq_size = 13)


library(ggplot2)

to_plot_df <- data.frame(Clique_size=c("1","==2",">=2",">=3",">=4",
                                       ">=5",">=6",">=7",">=8",">=9",
                                       ">=10",">=11",">=12",">=13"),
                percent=c(fr_orient_result_1_gre2$smallCliq_FR_percent,
                fr_orient_result_1_eq2$BigCliq_FR_percent,
                fr_orient_result_1_gre2$BigCliq_FR_percent,
                fr_orient_result_1_gre3$BigCliq_FR_percent,
                fr_orient_result_1_gre4$BigCliq_FR_percent,
                fr_orient_result_1_gre5$BigCliq_FR_percent,
                fr_orient_result_1_gre6$BigCliq_FR_percent,
                fr_orient_result_1_gre7$BigCliq_FR_percent,
                fr_orient_result_1_gre8$BigCliq_FR_percent,
                fr_orient_result_1_gre9$BigCliq_FR_percent,
                fr_orient_result_1_gre10$BigCliq_FR_percent,
                fr_orient_result_1_gre11$BigCliq_FR_percent,
                fr_orient_result_1_gre12$BigCliq_FR_percent,
                fr_orient_result_1_gre13$BigCliq_FR_percent)
                )

to_plot_df$Clique_size <- factor(to_plot_df$Clique_size, 
                                 levels = c("1","==2",">=2",">=3",">=4",
                                            ">=5",">=6",">=7",">=8",">=9",
                                            ">=10",">=11",">=12",">=13"))

ggplot(data = to_plot_df[c(1:11),], aes(x=Clique_size,y=percent)) + 
  geom_point(color="darkred",size=3) + 
  ylim(0,25) +
  geom_vline(xintercept = 3.5) +
  geom_hline(yintercept = (nrow(clic_orient[clic_orient$orient=="F,R",])/nrow(clic_orient)*100)) +
  theme_classic() + 
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 20)) 

ggsave("F3_IMR90FR_percent_vs_cliqueSize.pdf", plot = last_plot(), device = "pdf",
       scale = 1, width = 5, height = 5, units = "in",
       dpi = 300, limitsize = TRUE)

test_result_df <- data.frame(Clique_size=c("1","==2",">=2",">=3",">=4",
                                           ">=5",">=6",">=7",">=8",">=9",
                                           ">=10",">=11",">=12",">=13"),
                             percent=c(fr_orient_result_1_gre2$smallCliq_FR_percent,
                                       fr_orient_result_1_eq2$BigCliq_FR_percent,
                                       fr_orient_result_1_gre2$BigCliq_FR_percent,
                                       fr_orient_result_1_gre3$BigCliq_FR_percent,
                                       fr_orient_result_1_gre4$BigCliq_FR_percent,
                                       fr_orient_result_1_gre5$BigCliq_FR_percent,
                                       fr_orient_result_1_gre6$BigCliq_FR_percent,
                                       fr_orient_result_1_gre7$BigCliq_FR_percent,
                                       fr_orient_result_1_gre8$BigCliq_FR_percent,
                                       fr_orient_result_1_gre9$BigCliq_FR_percent,
                                       fr_orient_result_1_gre10$BigCliq_FR_percent,
                                       fr_orient_result_1_gre11$BigCliq_FR_percent,
                                       fr_orient_result_1_gre12$BigCliq_FR_percent,
                                       fr_orient_result_1_gre13$BigCliq_FR_percent),
                             pval=c("NA",fr_orient_result_1_eq2$statistics.p.value,
                                    fr_orient_result_1_gre2$statistics.p.value,
                                    fr_orient_result_1_gre3$statistics.p.value,
                                    fr_orient_result_1_gre4$statistics.p.value,
                                    fr_orient_result_1_gre5$statistics.p.value,
                                    fr_orient_result_1_gre6$statistics.p.value,
                                    fr_orient_result_1_gre7$statistics.p.value,
                                    fr_orient_result_1_gre8$statistics.p.value,
                                    fr_orient_result_1_gre9$statistics.p.value,
                                    fr_orient_result_1_gre10$statistics.p.value,
                                    fr_orient_result_1_gre11$statistics.p.value,
                                    fr_orient_result_1_gre12$statistics.p.value,
                                    fr_orient_result_1_gre13$statistics.p.value)
)

write.table(test_result_df, file = "f-test_stat_pval_and_percent_compared_with_equalto1.tsv",quote = F,
            sep = "\t",row.names = F, col.names = T)

test_result_df <- data.frame(Clique_size=c(">=2",">=3",">=4",
                                           ">=5",">=6",">=7",">=8",">=9",
                                           ">=10",">=11",">=12",">=13"),
                             
                             percent=c(fr_orient_result_2_gre2$BigCliq_FR_percent,
                                       fr_orient_result_2_gre3$BigCliq_FR_percent,
                                       fr_orient_result_2_gre4$BigCliq_FR_percent,
                                       fr_orient_result_2_gre5$BigCliq_FR_percent,
                                       fr_orient_result_2_gre6$BigCliq_FR_percent,
                                       fr_orient_result_2_gre7$BigCliq_FR_percent,
                                       fr_orient_result_2_gre8$BigCliq_FR_percent,
                                       fr_orient_result_2_gre9$BigCliq_FR_percent,
                                       fr_orient_result_2_gre10$BigCliq_FR_percent,
                                       fr_orient_result_2_gre11$BigCliq_FR_percent,
                                       fr_orient_result_2_gre12$BigCliq_FR_percent,
                                       fr_orient_result_2_gre13$BigCliq_FR_percent),
                             pval=c(fr_orient_result_2_gre2$statistics.p.value,
                                    fr_orient_result_2_gre3$statistics.p.value,
                                    fr_orient_result_2_gre4$statistics.p.value,
                                    fr_orient_result_2_gre5$statistics.p.value,
                                    fr_orient_result_2_gre6$statistics.p.value,
                                    fr_orient_result_2_gre7$statistics.p.value,
                                    fr_orient_result_2_gre8$statistics.p.value,
                                    fr_orient_result_2_gre9$statistics.p.value,
                                    fr_orient_result_2_gre10$statistics.p.value,
                                    fr_orient_result_2_gre11$statistics.p.value,
                                    fr_orient_result_2_gre12$statistics.p.value,
                                    fr_orient_result_2_gre13$statistics.p.value)
)

write.table(test_result_df, file = "f-test_stat_pval_and_percent_compared_with_equalto2.tsv",quote = F,
            sep = "\t",row.names = F, col.names = T)

fr_orient_result_1_eq2 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 1, big_cliq_size = 2,
                                              big_cliq_size_only_eq = TRUE)

#Total FR percent - no matter what category is that.
nrow(clic_orient[clic_orient$orient=="F,R",])/nrow(clic_orient)

# # of convergent tads
clic_orient_cp <- clic_orient
clic_orient_cp$sze_cat <- ifelse(clic_orient_cp$cliq_size>=8,">=8",clic_orient_cp$cliq_size)
freq_table <- as.data.frame(table(clic_orient_cp[clic_orient_cp$orient=="F,R","sze_cat"]))
freq_table$Var1 <- factor(freq_table$Var1,levels = c("1","2","3","4","5","6","7",">=8"))
ggplot(data = freq_table, aes(x=Var1,y=Freq)) + geom_bar(stat="identity", color="black", fill="white",size=1) + 
  theme_classic() +
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  xlab("Clique size") +
  ylab("Number of convergent TADs")

# # of tads
clic_orient_cp <- clic_orient
clic_orient_cp$sze_cat <- ifelse(clic_orient_cp$cliq_size>=8,">=8",clic_orient_cp$cliq_size)
freq_table <- as.data.frame(table(clic_orient_cp[,"sze_cat"]))
freq_table$Var1 <- factor(freq_table$Var1,levels = c("1","2","3","4","5","6","7",">=8"))
ggplot(data = freq_table, aes(x=Var1,y=Freq)) + geom_bar(stat="identity", color="black", fill="white",size=1) + 
  theme_classic() +
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  xlab("Clique size") +
  ylab("Number of all TADs")


fr_orient_result_1_gre2 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 1, big_cliq_size = 2)
fr_orient_result_1_eq2 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                             small_cliq_size = 1, big_cliq_size = 2,
                                             big_cliq_size_only_eq = TRUE)
fr_orient_result_1_gre3 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 1, big_cliq_size = 3)
fr_orient_result_1_gre4 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 1, big_cliq_size = 4)
fr_orient_result_1_gre5 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 1, big_cliq_size = 5)
fr_orient_result_1_gre6 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 1, big_cliq_size = 6)
fr_orient_result_1_gre7 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 1, big_cliq_size = 7)
fr_orient_result_1_gre8 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 1, big_cliq_size = 8)
fr_orient_result_1_gre9 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                              small_cliq_size = 1, big_cliq_size = 9)

fr_orient_result_1_gre10 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                               small_cliq_size = 1, big_cliq_size = 10)
fr_orient_result_1_gre11 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                               small_cliq_size = 1, big_cliq_size = 11)

fr_orient_result_1_gre12 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                               small_cliq_size = 1, big_cliq_size = 12)
fr_orient_result_1_gre13 <- run_FR_orient_test(clic_orient_df = clic_orient,
                                               small_cliq_size = 1, big_cliq_size = 13)


colnames(stats) <- c("cliq_cat","confident1","confident2","convergent_tads","total_tads_in_cat","global_convergent_tads","pvalue")
binom_test_for_cliq_against_global <- function(clic_orient)
{
  stats <- data.frame(cliq_cat=c(1:12),confident1=c(1:12),confident2=c(1:12),convergent_tads=c(1:12),
                      total_tads_in_cat=c(1:12),global_convergent_tads=c(1:12),pvalue=c(1:12))
  rownames(stats) <- c(">=2",">=3",">=4",">=5",">=6",">=7",">=8",">=9",">=10",">=11",">=12",">=13")
  global_avg <- nrow(clic_orient[clic_orient$orient=="F,R",])/nrow(clic_orient)
  for (i in c(2:13)) 
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
