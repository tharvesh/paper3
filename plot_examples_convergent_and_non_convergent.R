
m <-read.table("chr1_Ex1_plot_cliq_12.mat",colClasses = "numeric")
rownames(m) <- seq(207800000,207800000+(nrow(m)*5000)-5000,by=5000)
colnames(m) <- seq(207800000,207800000+(nrow(m)*5000)-5000,by=5000)
m <- as.matrix(m)

image(as.numeric(rownames(m)),as.numeric(rownames(m)),t(log(m)),col = colorRampPalette(c("white", "red"))(25),las = 1,xlab = "",ylab = "",xaxt="n",yaxt = "n")
axis(1, at=seq(207800000,209750000,by=100000),labels = seq(207800000,209750000,by=100000)/1000000, las=1, cex.axis=2)

m2 <- read.table("chr18_Ex2_plot_cliq_12.mat",colClasses = "numeric")
rownames(m2) <- seq(70050000,70050000+(nrow(m2)*5000)-5000,by=5000)
colnames(m2) <- seq(70050000,70050000+(nrow(m2)*5000)-5000,by=5000)
m2 <- as.matrix(m2)
image(as.numeric(rownames(m2)),as.numeric(rownames(m2)),t(log(m2)),col = colorRampPalette(c("white", "red"))(25),las = 1,xlab = "",ylab = "",xaxt="n",yaxt = "n")
axis(1, at=seq(70050000,72800000,by=100000),labels = seq(70050000,72800000,by=100000)/1000000, las=1, cex.axis=2)

m3 <- read.table("chr14_Ex3_plot_cliq_12.mat",colClasses = "numeric")
rownames(m3) <- seq(39000000,39000000+(nrow(m3)*5000)-5000,by=5000)
colnames(m3) <- seq(39000000,39000000+(nrow(m3)*5000)-5000,by=5000)
m3 <- as.matrix(m3)
image(as.numeric(rownames(m3)),as.numeric(rownames(m3)),t(log(m3)),col = colorRampPalette(c("white", "red"))(25),las = 1,xlab = "",ylab = "",xaxt="n",yaxt = "n")
axis(1, at=seq(39000000,42600000,by=100000),labels = seq(39000000,42600000,by=100000)/1000000, las=1, cex.axis=2)
