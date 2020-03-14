library(Matrix)

rescale <- function(x, newrange=range(x)){
  xrange <- range(x)
  mfac <- (newrange[2]-newrange[1])/(xrange[2]-xrange[1])
  newrange[1]+(x-xrange[1])*mfac
}

ResizeMat <- function(mat, ndim=dim(mat)){
  if(!require(fields)) stop("`fields` required.")
  
  # input object
  odim <- dim(mat)
  obj <- list(x= 1:odim[1], y=1:odim[2], z= mat)
  
  # output object
  ans <- matrix(NA, nrow=ndim[1], ncol=ndim[2])
  ndim <- dim(ans)
  
  # rescaling
  ncord <- as.matrix(expand.grid(seq_len(ndim[1]), seq_len(ndim[2])))
  loc <- ncord
  loc[,1] = rescale(ncord[,1], c(1,odim[1]))
  loc[,2] = rescale(ncord[,2], c(1,odim[2]))
  
  # interpolation
  ans[ncord] <- interp.surface(obj, loc)
  
  ans
}

library(OpenImageR)

aggregated_mean_of_tads <- function(path_to_mat_files,target_size=25) 
{
  require(OpenImageR)
  filenames <- list.files(path_to_mat_files,full.names = T)
  lv <- list()
  i=1
  
  for (f in filenames) 
  {
    m <- as.matrix(read.table(f,colClasses = "numeric"))
    g <- m/sum(m)
    if (dim(m)[1]==dim(m)[2]) 
    {
      if (sum(diag(m))==0) 
      {
        if(sum(m)!=0)
        {
          res <- resizeImage(g,target_size,target_size,"nearest")
          lv[[i]] <- res
        }
      }
      i=i+1
    }
  }
  
  if(length(unique(sapply(lv, is.null)))>1)
  {
    lv2 <- lv[-which(sapply(lv, is.null))]
  }
  else
  {
    lv2 <- lv
  }
  
  removing_ind <- c()
  ind=0
  for (mat in lv2) 
  {
    ind=ind+1
    if(dim(mat)[1]!=target_size | dim(mat)[2]!=target_size )
    {
      removing_ind <- c(removing_ind,ind)
    }
  }
  
  if(length(removing_ind)>0)
  {
    lv2 <- lv2[-removing_ind]
  }
  
  mean = Reduce("+", lv2) / length(lv2)
  return(mean)
}

aggregated_mean_of_tads_size_selected <- function(path_to_mat_files,target_size=25,min_size_tad_bins=50,max_size_tad_bins=300) 
{
  require(OpenImageR)
  filenames <- list.files(path_to_mat_files,full.names = T)
  lv <- list()
  i=1
  
  for (f in filenames) 
  {
    m <- as.matrix(read.table(f,colClasses = "numeric"))
    g <- m/sum(m)
    if (dim(m)[1]==dim(m)[2]) 
    {
      if (sum(diag(m))==0) 
      {
        if(sum(m)!=0)
        {
          if (nrow(g)>=min_size_tad_bins & nrow(g) <= max_size_tad_bins) {
            res <- resizeImage(g,target_size,target_size,"nearest")
            lv[[i]] <- res 
          }
        }
      }
      i=i+1
    }
  }
  
  if(length(unique(sapply(lv, is.null)))>1)
  {
    lv2 <- lv[-which(sapply(lv, is.null))]
  }
  else
  {
    lv2 <- lv
  }
  
  removing_ind <- c()
  ind=0
  for (mat in lv2) 
  {
    ind=ind+1
    if(dim(mat)[1]!=target_size | dim(mat)[2]!=target_size )
    {
      removing_ind <- c(removing_ind,ind)
    }
  }
  
  if(length(removing_ind)>0)
  {
    lv2 <- lv2[-removing_ind]
  }
  
  mean = Reduce("+", lv2) / length(lv2)
  return(mean)
}

write_mean_of_reduced_matrices <- function(path_big_cliq,target_size=25)
{
  mean_mat <- aggregated_mean_of_tads(path_big_cliq,target_size = target_size)
  out_path <- paste(basename(path_big_cliq),"_mean_of_reduced_mat.tsv",sep = "")
  write.table(mean_mat,out_path,col.names = F,sep = "\t",row.names = F,quote = F)
}


write_log_ratio_bw_mean_tads <- function(path_big_cliq,path_cliq_size_small_file,target_size=25,log_max=8)
{
  require(tools)
  mean_small <- as.matrix(read.table(path_cliq_size_small_file,colClasses = "numeric"))
  mean_big <- aggregated_mean_of_tads(path_big_cliq,target_size = target_size)
  plt <- log2(mean_small) - log2(mean_big)
  out_path <- paste("log2",basename(path_big_cliq),"_diff_",file_path_sans_ext(path_cliq_size_small_file),".tsv",sep = "")
  write.table(plt,out_path,col.names = F,sep = "\t",row.names = F,quote = F)
  #levelplot(plt,at=seq(0,8,length=20),col.regions=coolwarm(1000),region=T,xlab="",ylab="",column.values = seq_len(nrow(plt)),scales = list(tck = c(0,0)),
  #          labels=NULL)
}

aggregated_mean_of_tads_no_norm <- function(path_to_mat_files,target_size=25) 
{
  require(OpenImageR)
  filenames <- list.files(path_to_mat_files,full.names = T)
  lv <- list()
  i=1
  
  for (f in filenames) 
  {
    m <- as.matrix(read.table(f,colClasses = "numeric"))
    if (dim(m)[1]==dim(m)[2]) 
    {
      if (sum(diag(m))==0) 
      {
        if(sum(m)!=0)
        {
          res <- resizeImage(m,target_size,target_size,"nearest")
          lv[[i]] <- res
        }
      }
      i=i+1
    }
  }
  
  if(length(unique(sapply(lv, is.null)))>1)
  {
    lv2 <- lv[-which(sapply(lv, is.null))]
  }
  else
  {
    lv2 <- lv
  }
  
  removing_ind <- c()
  ind=0
  for (mat in lv2) 
  {
    ind=ind+1
    if(dim(mat)[1]!=target_size | dim(mat)[2]!=target_size )
    {
      removing_ind <- c(removing_ind,ind)
    }
  }
  
  if(length(removing_ind)>0)
  {
    lv2 <- lv2[-removing_ind]
  }
  
  mean = Reduce("+", lv2) / length(lv2)
  return(mean)
}



plot_log_ratio_bw_mean_tads <- function(path_big_cliq,path_small_cliq,target_size=25)
{
  require(rafalib)
  mean_small <- aggregated_mean_of_tads_size_selected(path_small_cliq,target_size = target_size)
  mean_big <- aggregated_mean_of_tads_size_selected(path_big_cliq,target_size = target_size)
  plt <- log2(mean_small) - log2(mean_big)
  #pdf_file <- paste(basename(path_big_cliq),"_diff_",file_path_sans_ext(path_cliq_size_small_file),".pdf",sep = "")
  #levelplot(plt,cut=20,col.regions=coolwarm(1000),region=T,xlab="",ylab="",column.values = seq_len(nrow(plt)),scales = list(tck = c(0,0)),
  #          labels=NULL)
  mypar(1,3)
  imagemat(log2(mean_small), col = colorRampPalette(c("white", "red"))(25),main=basename(path_small_cliq))
  imagemat(log2(mean_big), col = colorRampPalette(c("white", "red"))(25),main=basename(path_big_cliq))
  imagemat(plt, col = colorRampPalette(c("white", "red"))(25),main="ratio")
}

plot_mean_tads <- function(path_cliq_folder,main,target_size=25,zlim=c(0,0.3))
{
  require(rafalib)
  mean <- aggregated_mean_of_tads_no_norm(path_cliq_folder,target_size = target_size)
  image(mean, col = colorRampPalette(c("white", "red"))(25),main=main,las = 1,xlab = "",ylab = "",xaxt="n",yaxt = "n",zlim = zlim)
}

plot_log_ratio_bw_mean_tads_old <- function(path_big_cliq,path_cliq_size_small_file,target_size=25)
{
  require(tools)
  mean_small <- as.matrix(read.table(path_cliq_size_small_file,colClasses = "numeric"))
  mean_big <- aggregated_mean_of_tads(path_big_cliq,target_size = target_size)
  plt <- log2(mean_small) - log2(mean_big)
  pdf_file <- paste(basename(path_big_cliq),"_diff_",file_path_sans_ext(path_cliq_size_small_file),".pdf",sep = "")
  levelplot(plt,cut=20,col.regions=coolwarm(1000),region=T,xlab="",ylab="",column.values = seq_len(nrow(plt)),scales = list(tck = c(0,0)),
            labels=NULL)
}

plot_log_ratio_bw_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq2",path_cliq_size_small_file = "cliq_size_1.tsv",
                            target_size = 25)

plot_log_ratio_bw_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq3",path_cliq_size_small_file = "cliq_size_1.tsv",
                            target_size = 25)
plot_log_ratio_bw_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq4",path_cliq_size_small_file = "cliq_size_1.tsv",
                            target_size = 25)
plot_log_ratio_bw_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq5",path_cliq_size_small_file = "cliq_size_1.tsv",
                            target_size = 25)
plot_log_ratio_bw_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq6",path_cliq_size_small_file = "cliq_size_1.tsv",
                            target_size = 25)
plot_log_ratio_bw_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq7",path_cliq_size_small_file = "cliq_size_1.tsv",
                            target_size = 25)
plot_log_ratio_bw_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq8",path_cliq_size_small_file = "cliq_size_1.tsv",
                            target_size = 25)
plot_log_ratio_bw_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq9",path_cliq_size_small_file = "cliq_size_1.tsv",
                            target_size = 25)
plot_log_ratio_bw_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq10",path_cliq_size_small_file = "cliq_size_1.tsv",
                            target_size = 25)

write_log_ratio_bw_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq2",path_cliq_size_small_file = "cliq_size_1.tsv",
                            target_size = 25)

write_log_ratio_bw_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq3",path_cliq_size_small_file = "cliq_size_1.tsv",
                            target_size = 25)
write_log_ratio_bw_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq4",path_cliq_size_small_file = "cliq_size_1.tsv",
                            target_size = 25)
write_log_ratio_bw_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq5",path_cliq_size_small_file = "cliq_size_1.tsv",
                            target_size = 25)
write_log_ratio_bw_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq6",path_cliq_size_small_file = "cliq_size_1.tsv",
                            target_size = 25)
write_log_ratio_bw_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq7",path_cliq_size_small_file = "cliq_size_1.tsv",
                            target_size = 25)
write_log_ratio_bw_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq8",path_cliq_size_small_file = "cliq_size_1.tsv",
                            target_size = 25)
write_log_ratio_bw_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq9",path_cliq_size_small_file = "cliq_size_1.tsv",
                            target_size = 25)
write_log_ratio_bw_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq10",path_cliq_size_small_file = "cliq_size_1.tsv",
                            target_size = 25)


write_mean_of_reduced_matrices("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq2")
write_mean_of_reduced_matrices("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq3")
write_mean_of_reduced_matrices("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq4")
write_mean_of_reduced_matrices("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq5")
write_mean_of_reduced_matrices("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq6")
write_mean_of_reduced_matrices("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq7")
write_mean_of_reduced_matrices("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq8")
write_mean_of_reduced_matrices("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq9")
write_mean_of_reduced_matrices("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq10")



path <- "/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/"
for (ind in c(1:12,14)) 
{
  folder <- paste("cliq_size_",ind,sep="")
  folderpath <- paste(path,folder,sep="")
  outpath <- paste(folder,".tsv",sep = "")
  cat(outpath,"\n")
  #mean <- aggregated_mean_of_tads(folderpath,target_size = 25)
  #write.table(mean,outpath,col.names = F,sep = "\t",row.names = F,quote = F)
}


par(mar=c(1, 1, 1, 1))
par(mfrow=c(3,3))
plot_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_1",main="Clique size = 1")
plot_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_2",main="Clique size = 2")
plot_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq2",main="Clique size >= 2")
plot_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq3",main="Clique size >= 3")
plot_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq4",main="Clique size >= 4")
plot_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq5",main="Clique size >= 5")
plot_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq6",main="Clique size >= 6")
plot_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq7",main="Clique size >= 7")
plot_mean_tads("/Users/tmali/XX_cliques/CTCF_motif_analysis/triangle_hist/matrix_format/cliq_size_geq8",main="Clique size >= 8")
