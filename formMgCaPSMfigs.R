#Figure for foram Mg/Ca-environmental data calibration

setwd("/Volumes/GoogleDrive/My Drive/R_documents/PSR");

library(ncdf4)
library(maps);library(mapproj); library(mapdata);

#coretop foram MgCa data
MgCa.dat.Npachy<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/ForamMg_Npachys.csv")
MgCa.dat.Gruber<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/ForamMg_Gruber.csv")
MgCa.dat.Ginflata<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/ForamMg_Ginflata.csv")
MgCa.dat.Gbulloides<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/ForamMg_Gbulloides.csv")

#results from iterative regressions
CE.Npachy<-read.csv("Npachy_summary")
CE.Gruber<-read.csv("Gruber_summary")
CE.Ginflata<-read.csv("Ginflata_summary")
CE.Gbulloides<-read.csv("Gbulloides_summary")
CE.all<-read.csv("all_summary")

#Figure 1 map of spatial distribution 
dev.new(width=12, height=6);
world<-map('worldHires',interior=FALSE,ylim=c(-88,88),xlim=c(-180,180))
points(MgCa.dat.Gruber$lon,MgCa.dat.Gruber$lat,cex=0.8,col="cornflowerblue")
points(MgCa.dat.Ginflata$lon,MgCa.dat.Ginflata$lat,cex=1.1,col="firebrick")
points(MgCa.dat.Gbulloides$lon,MgCa.dat.Gbulloides$lat,cex=1.4,col="gold3")
points(MgCa.dat.Npachy$lon,MgCa.dat.Npachy$lat,cex=0.5,col="palegreen3")
legend(30,88,legend=c("N. pachyderma","G. ruber","G. inflata","G. bulloides"), pch=c(1,1,1,1), col=c("palegreen3","cornflowerblue","firebrick","gold3"),bg="white")
dev.print(pdf,"MgCa_map.pdf");
dev.off();

#Figure 2 summary of CE
dev.new(width=8, height=11);
par(mfrow=c(2,2))
par(mar=c(15,5,2,2))
plot(CE.Npachy$CE,xaxt="n",pch=16,xlab="",las=1,ylim=c(-0.1,1),ylab="CE",main="N. pachyderma")
axis(1,at=1:nrow(CE.Npachy),labels=CE.Npachy$variables,las=2,cex=0.5)
segments(seq(1:nrow(CE.Npachy)),CE.Npachy$CE-CE.Npachy$CE_95ci,seq(1:nrow(CE.Npachy)),CE.Npachy$CE+CE.Npachy$CE_95ci)
abline(0,0)
plot(CE.Gruber$CE,xaxt="n",pch=16,xlab="",las=1,ylim=c(-0.1,1),ylab="CE",main="G. ruber")
axis(1,at=1:nrow(CE.Gruber),labels=CE.Gruber$variables,las=2,cex=0.5)
segments(seq(1:nrow(CE.Gruber)),CE.Gruber$CE-CE.Gruber$CE_95ci,seq(1:nrow(CE.Gruber)),CE.Gruber$CE+CE.Gruber$CE_95ci)
abline(0,0)
plot(CE.Ginflata$CE,xaxt="n",pch=16,xlab="",las=1,ylim=c(-0.1,1),ylab="CE",main="G. inflata")
axis(1,at=1:nrow(CE.Ginflata),labels=CE.Ginflata$variables,las=2,cex=0.5)
segments(seq(1:nrow(CE.Ginflata)),CE.Ginflata$CE-CE.Ginflata$CE_95ci,seq(1:nrow(CE.Ginflata)),CE.Ginflata$CE+CE.Ginflata$CE_95ci)
abline(0,0)
plot(CE.Gbulloides$CE,xaxt="n",pch=16,xlab="",las=1,ylim=c(-0.1,1),ylab="CE",main="G. bulloides")
axis(1,at=1:nrow(CE.Gbulloides),labels=CE.Gbulloides$variables,las=2,cex=0.5)
segments(seq(1:nrow(CE.Gbulloides)),CE.Gbulloides$CE-CE.Gbulloides$CE_95ci,seq(1:nrow(CE.Gbulloides)),CE.Gbulloides$CE+CE.Gbulloides$CE_95ci)
abline(0,0)

dev.copy(pdf,"CEplot_seasons.pdf")
dev.off()

#Figure 3 summary of RMSE
dev.new(width=11, height=8);
par(mfrow=c(2,3))
par(mar=c(15,5,2,2))
plot(CE.Npachy$RMSE,xaxt="n",pch=16,xlab="",las=1,ylim=c(0,1.4),ylab="RMSE",main="N. pachyderma")
axis(1,at=1:nrow(CE.Npachy),labels=CE.Npachy$variables,las=2,cex=0.5)
segments(seq(1:nrow(CE.Npachy)),CE.Npachy$RMSE-CE.Npachy$RMSE_95ci,seq(1:nrow(CE.Npachy)),CE.Npachy$RMSE+CE.Npachy$RMSE_95ci)
abline(0,0)
points(winterCE.Npachy$RMSE,xaxt="n",col="cornflowerblue")
points(springCE.Npachy$RMSE,xaxt="n",col="palegreen3")
points(summerCE.Npachy$RMSE,xaxt="n",col="firebrick")
points(fallCE.Npachy$RMSE,xaxt="n",col="gold3")
legend("topleft",legend=c("annual","winter","spring","summer","fall"), pch=c(16,1,1,1,1), col=c("black","cornflowerblue","palegreen3","firebrick","gold3"),bty="n")

plot(CE.Gruber$RMSE,xaxt="n",pch=16,xlab="",las=1,ylim=c(0,1.4),ylab="RMSE",main="G. ruber")
axis(1,at=1:nrow(CE.Gruber),labels=CE.Gruber$variables,las=2,cex=0.5)
segments(seq(1:nrow(CE.Gruber)),CE.Gruber$RMSE-CE.Gruber$RMSE_95ci,seq(1:nrow(CE.Gruber)),CE.Gruber$RMSE+CE.Gruber$RMSE_95ci)
abline(0,0)
points(winterCE.Gruber$RMSE,xaxt="n",col="cornflowerblue")
points(springCE.Gruber$RMSE,xaxt="n",col="palegreen3")
points(summerCE.Gruber$RMSE,xaxt="n",col="firebrick")
points(fallCE.Gruber$RMSE,xaxt="n",col="gold3")

plot(CE.Ginflata$RMSE,xaxt="n",pch=16,xlab="",las=1,ylim=c(0,1.4),ylab="RMSE",main="G. inflata")
axis(1,at=1:nrow(CE.Ginflata),labels=CE.Ginflata$variables,las=2,cex=0.5)
segments(seq(1:nrow(CE.Ginflata)),CE.Ginflata$RMSE-CE.Ginflata$RMSE_95ci,seq(1:nrow(CE.Ginflata)),CE.Ginflata$RMSE+CE.Ginflata$RMSE_95ci)
abline(0,0)
points(winterCE.Ginflata$RMSE,xaxt="n",col="cornflowerblue")
points(springCE.Ginflata$RMSE,xaxt="n",col="palegreen3")
points(summerCE.Ginflata$RMSE,xaxt="n",col="firebrick")
points(fallCE.Ginflata$RMSE,xaxt="n",col="gold3")

plot(CE.Gbulloides$RMSE,xaxt="n",pch=16,xlab="",las=1,ylim=c(0,2.5),ylab="RMSE",main="G. bulloides")
axis(1,at=1:nrow(CE.Gbulloides),labels=CE.Gbulloides$variables,las=2,cex=0.5)
segments(seq(1:nrow(CE.Gbulloides)),CE.Gbulloides$RMSE-CE.Gbulloides$RMSE_95ci,seq(1:nrow(CE.Gbulloides)),CE.Gbulloides$RMSE+CE.Gbulloides$RMSE_95ci)
abline(0,0)
points(winterCE.Gbulloides$RMSE,xaxt="n",col="cornflowerblue")
points(springCE.Gbulloides$RMSE,xaxt="n",col="palegreen3")
points(summerCE.Gbulloides$RMSE,xaxt="n",col="firebrick")
points(fallCE.Gbulloides$RMSE,xaxt="n",col="gold3")

plot(CE.all$RMSE,xaxt="n",pch=16,xlab="",las=1,ylim=c(0,1.4),ylab="RMSE",main="all species")
axis(1,at=1:nrow(CE.all),labels=CE.all$variables,las=2,cex=0.5)
segments(seq(1:nrow(CE.all)),CE.all$RMSE-CE.all$RMSE_95ci,seq(1:nrow(CE.all)),CE.all$RMSE+CE.all$RMSE_95ci)
abline(0,0)
points(winterCE.all$RMSE,xaxt="n",col="cornflowerblue")
points(springCE.all$RMSE,xaxt="n",col="palegreen3")
points(summerCE.all$RMSE,xaxt="n",col="firebrick")
points(fallCE.all$RMSE,xaxt="n",col="gold3")

dev.copy(pdf,"RMSEplot_seasons.pdf")
dev.off()

#Figure 4 summary of r2
dev.new(width=11, height=8);
par(mfrow=c(2,3))
par(mar=c(15,5,2,2))
plot(CE.Npachy$r2,xaxt="n",pch=16,xlab="",las=1,ylim=c(0,1),ylab="r2",main="N. pachyderma")
axis(1,at=1:nrow(CE.Npachy),labels=CE.Npachy$variables,las=2,cex=0.5)
segments(seq(1:nrow(CE.Npachy)),CE.Npachy$r2-CE.Npachy$r2_95ci,seq(1:nrow(CE.Npachy)),CE.Npachy$r2+CE.Npachy$r2_95ci)
abline(0,0)
points(winterCE.Npachy$r2,xaxt="n",col="cornflowerblue")
points(springCE.Npachy$r2,xaxt="n",col="palegreen3")
points(summerCE.Npachy$r2,xaxt="n",col="firebrick")
points(fallCE.Npachy$r2,xaxt="n",col="gold3")
legend("topleft",legend=c("annual","winter","spring","summer","fall"), pch=c(16,1,1,1,1), col=c("black","cornflowerblue","palegreen3","firebrick","gold3"),bty="n")

plot(CE.Gruber$r2,xaxt="n",pch=16,xlab="",las=1,ylim=c(0,1),ylab="r2",main="G. ruber")
axis(1,at=1:nrow(CE.Gruber),labels=CE.Gruber$variables,las=2,cex=0.5)
segments(seq(1:nrow(CE.Gruber)),CE.Gruber$r2-CE.Gruber$r2_95ci,seq(1:nrow(CE.Gruber)),CE.Gruber$r2+CE.Gruber$r2_95ci)
abline(0,0)
points(winterCE.Gruber$r2,xaxt="n",col="cornflowerblue")
points(springCE.Gruber$r2,xaxt="n",col="palegreen3")
points(summerCE.Gruber$r2,xaxt="n",col="firebrick")
points(fallCE.Gruber$r2,xaxt="n",col="gold3")

plot(CE.Ginflata$r2,xaxt="n",pch=16,xlab="",las=1,ylim=c(0,1),ylab="r2",main="G. inflata")
axis(1,at=1:nrow(CE.Ginflata),labels=CE.Ginflata$variables,las=2,cex=0.5)
segments(seq(1:nrow(CE.Ginflata)),CE.Ginflata$r2-CE.Ginflata$r2_95ci,seq(1:nrow(CE.Ginflata)),CE.Ginflata$r2+CE.Ginflata$r2_95ci)
abline(0,0)
points(winterCE.Ginflata$r2,xaxt="n",col="cornflowerblue")
points(springCE.Ginflata$r2,xaxt="n",col="palegreen3")
points(summerCE.Ginflata$r2,xaxt="n",col="firebrick")
points(fallCE.Ginflata$r2,xaxt="n",col="gold3")

plot(CE.Gbulloides$r2,xaxt="n",pch=16,xlab="",las=1,ylim=c(0,1),ylab="r2",main="G. bulloides")
axis(1,at=1:nrow(CE.Gbulloides),labels=CE.Gbulloides$variables,las=2,cex=0.5)
segments(seq(1:nrow(CE.Gbulloides)),CE.Gbulloides$r2-CE.Gbulloides$r2_95ci,seq(1:nrow(CE.Gbulloides)),CE.Gbulloides$r2+CE.Gbulloides$r2_95ci)
abline(0,0)
points(winterCE.Gbulloides$r2,xaxt="n",col="cornflowerblue")
points(springCE.Gbulloides$r2,xaxt="n",col="palegreen3")
points(summerCE.Gbulloides$r2,xaxt="n",col="firebrick")
points(fallCE.Gbulloides$r2,xaxt="n",col="gold3")

plot(CE.all$r2,xaxt="n",pch=16,xlab="",las=1,ylim=c(0,1),ylab="r2",main="all species")
axis(1,at=1:nrow(CE.all),labels=CE.all$variables,las=2,cex=0.5)
segments(seq(1:nrow(CE.all)),CE.all$r2-CE.all$r2_95ci,seq(1:nrow(CE.all)),CE.all$r2+CE.all$r2_95ci)
abline(0,0)
points(winterCE.all$r2,xaxt="n",col="cornflowerblue")
points(springCE.all$r2,xaxt="n",col="palegreen3")
points(summerCE.all$r2,xaxt="n",col="firebrick")
points(fallCE.all$r2,xaxt="n",col="gold3")

dev.copy(pdf,"r2plot_seasons.pdf")
dev.off()

#Figure X comparison with originally published SST
DT<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/O2K_DT.csv")
DT.stda.bin<-matrix(,nrow=1,ncol=4)
dev.new(width=8, height=5);
par(mfrow=c(1,2))
plot(DT$temp[which(DT$Species=="G.ruber")],DT$Tonly_cps[which(DT$Species=="G.ruber")],xlab="T (ºC, original)",ylab="T (ºC, this study)", xlim=c(0,32),ylim=c(0,32),las=1,col="cornflowerblue",main="T-only",pch=20,cex=0.4)
points(DT$temp[which(DT$Species=="G.bulloides")],DT$Tonly_cps[which(DT$Species=="G.bulloides")], col="gold3",pch=20,cex=0.4)
points(DT$temp[which(DT$Species=="G.inflata")],DT$Tonly_cps[which(DT$Species=="G.inflata")], col="firebrick",pch=20,cex=0.4)
points(DT$temp[which(DT$Species=="N.pachyderma")],DT$Tonly_cps[which(DT$Species=="N.pachyderma")], col="palegreen3",pch=20,cex=0.4)
abline(0,1)
pal<-c("gold3","firebrick","cornflowerblue","palegreen3")


plot(DT$temp[which(DT$Species=="G.ruber")],DT$Tmult_cps[which(DT$Species=="G.ruber")],xlab="T (ºC, original)",ylab="T (ºC, this study)", xlim=c(0,32),ylim=c(0,32),las=1,col="cornflowerblue",main="multivariate",pch=20,cex=0.4)
points(DT$temp[which(DT$Species=="G.bulloides")],DT$Tmult_cps[which(DT$Species=="G.bulloides")], col="gold3",pch=20,cex=0.4)
points(DT$temp[which(DT$Species=="G.inflata")],DT$Tmult_cps[which(DT$Species=="G.inflata")], col="firebrick",pch=20,cex=0.4)
points(DT$temp[which(DT$Species=="N.pachyderma")],DT$Tmult_cps[which(DT$Species=="N.pachyderma")], col="palegreen3",pch=20,cex=0.4)
abline(0,1)
pal<-c("gold3","firebrick","cornflowerblue","palegreen3")
legend("topleft",legend=levels(DT$Species),text.col=pal[seq(1,4)],bty="n",cex=0.8)


dev.new(width=6, height=6);
plot(0,type="n",xlab="mean change in T (ºC, this study-original)",ylab="change in std. dev. (ºC, this study-original)",xlim=c(-5,5),ylim=c(-1,1.5),las=1)
abline(h=0,col="grey50");abline(v=0,col="grey50")
for (i in 1:length(unique(DT$index))) {
	tmp<-DT[which(DT$index==unique(DT$index)[i]),]
	var<-sd(tmp$temp_cps)-sd(tmp$temp)
	Dtemp<-mean(tmp$temp_cps-tmp$temp)
	text(Dtemp,var,labels=tmp$index,col=pal[as.numeric(tmp$Species[1])])
	temp.sa<-(tmp$temp-mean(tmp$temp))/sd(tmp$temp)											#original T standard anomaly
	temp_cps.sa<-(tmp$temp_cps-mean(tmp$temp_cps))/sd(tmp$temp_cps)							#this study T standard anomaly					
	temp.sa.bin<-as.numeric(tapply(temp.sa,cut(tmp$year,seq(0,2000,by=200)),mean))			#bin every 200 years
	temp_cps.sa.bin<-as.numeric(tapply(temp_cps.sa,cut(tmp$year,seq(0,2000,by=200)),mean))
	tmp1<-cbind(rep(tmp$index[1],10),seq(100,1900,by=200),temp.sa.bin,temp_cps.sa.bin)
	DT.stda.bin<-rbind(DT.stda.bin,tmp1)
	}
legend("bottomleft",legend=levels(DT$Species),text.col=pal[seq(1,4)],bty="n")
DT.stda.bin<-DT.stda.bin[-1,]
dev.new(width=7, height=6);
par(mfrow=c(2,1),mar=c(5,5,1,1))
boxplot(DT.stda.bin[,3]~DT.stda.bin[,2],names=seq(100,1900,200),las=1,xlab="",ylab="standardized T anomaly",ylim=c(-2.5,2), main="original")
boxplot(DT.stda.bin[,4]~DT.stda.bin[,2],names=seq(100,1900,200),las=1,xlab="center year of bin (A.D.)",ylab="standardized T anomaly",ylim=c(-2.5,2),main="this study")
