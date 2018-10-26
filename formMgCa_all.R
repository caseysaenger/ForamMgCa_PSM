#calibrate and validate forward models of foraminifera MgCa using an interative approach with CE as metric 

setwd("/Volumes/GoogleDrive/My Drive/R_documents/PSR");

#Resutls from formMgCa_cal.val suggest the following models: 
#N. pachyderma: f(T,size)
#G. ruber: f(T), f(T,size,OmegaCdeep) 
#G. inflata: f(T,pHsurf)
#G. bulloides: f(T) f(T,OmegaCdeep,size)
#Construct models using the variables above and ALL Mg/Ca for a species

source("/Volumes/GoogleDrive/My Drive/R_documents/funcs/cal.val.R")
source("/Volumes/GoogleDrive/My Drive/R_documents/funcs/cal.val.nls.R")

library(ncdf4)
library(leaps)
library(abind)
library(plyr)

#constants
OmegaCcut<-1.3																				#threshold below which calcite saturation has an effect

#coretop foram MgCa data
MgCa.dat.Npachy<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/ForamMg_Npachys.csv")
MgCa.dat.Gruber<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/ForamMg_Gruber.csv")
MgCa.dat.Ginflata<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/ForamMg_Ginflata.csv")
MgCa.dat.Gbulloides<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/ForamMg_Gbulloides.csv")

#culture data
MgCa.dat.culture<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/ForamMg_culture.csv")

#empirically account for reductive cleaning following Rosenthal et al., 2004 relationship
red<-which(MgCa.dat.Npachy$reductive=="Y")
MgCa.dat.Npachy$Mg.Ca[red]<-MgCa.dat.Npachy$Mg.Ca[red]/0.85									#omit intercept for N. pachy as this brings Yu, 2008 and E+G, 2000 into better agreement
red<-which(MgCa.dat.Gruber$reductive=="Y")
MgCa.dat.Gruber$Mg.Ca[red]<-(MgCa.dat.Gruber$Mg.Ca[red]-0.2)/0.85
red<-which(MgCa.dat.Ginflata$reductive=="Y")
MgCa.dat.Ginflata$Mg.Ca[red]<-(MgCa.dat.Ginflata$Mg.Ca[red]-0.2)/0.85
red<-which(MgCa.dat.Gbulloides$reductive=="Y")
MgCa.dat.Gbulloides$Mg.Ca[red]<-(MgCa.dat.Gbulloides$Mg.Ca[red]-0.2)/0.85

#envrionmental data
#salinity
nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/WOA/woa13_decav_s00_01v2.nc")
s.ann<-ncvar_get(nc,"s_an",start=c(1,1,1,1),count=c(-1,-1,40,-1))	
depth<-nc$dim$depth$vals[1:40]
lat<-nc$dim$lat$vals
lon<-nc$dim$lon$vals
rm(nc)
nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/WOA/woa13_decav_s13_01v2.nc")
s.winter.tmp<-ncvar_get(nc,"s_an",start=c(1,1,1,1),count=c(-1,-1,40,-1))
s.winterNH<-s.winter.tmp[,which(lat>=0),]
s.summerSH<-s.winter.tmp[,which(lat<0),]
rm(nc)
nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/WOA/woa13_decav_s15_01v2.nc")
s.summer.tmp<-ncvar_get(nc,"s_an",start=c(1,1,1,1),count=c(-1,-1,40,-1))
s.summerNH<-s.summer.tmp[,which(lat>=0),]
s.winterSH<-s.summer.tmp[,which(lat<0),]
rm(nc)
s.winter<-abind(s.winterSH,s.winterNH,along=2)
s.summer<-abind(s.summerSH,s.summerNH,along=2)

nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/WOA/woa13_decav_s14_01v2.nc")
s.spring.tmp<-ncvar_get(nc,"s_an",start=c(1,1,1,1),count=c(-1,-1,40,-1))
s.springNH<-s.spring.tmp[,which(lat>=0),]
s.fallSH<-s.spring.tmp[,which(lat<0),]
rm(nc)
nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/WOA/woa13_decav_s16_01v2.nc")
s.fall.tmp<-ncvar_get(nc,"s_an",start=c(1,1,1,1),count=c(-1,-1,40,-1))
s.fallNH<-s.fall.tmp[,which(lat>=0),]
s.springSH<-s.fall.tmp[,which(lat<0),]
rm(nc)
s.spring<-abind(s.springSH,s.springNH,along=2)
s.fall<-abind(s.fallSH,s.fallNH,along=2)

#temperature
nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/WOA/woa13_decav_t00_01v2.nc")
t.ann<-ncvar_get(nc,"t_an",start=c(1,1,1,1),count=c(-1,-1,40,-1))	
rm(nc)
nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/WOA/woa13_decav_t13_01v2.nc")
t.winter.tmp<-ncvar_get(nc,"t_an",start=c(1,1,1,1),count=c(-1,-1,40,-1))
t.winterNH<-t.winter.tmp[,which(lat>=0),]
t.summerSH<-t.winter.tmp[,which(lat<0),]
rm(nc)
nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/WOA/woa13_decav_t15_01v2.nc")
t.summer.tmp<-ncvar_get(nc,"t_an",start=c(1,1,1,1),count=c(-1,-1,40,-1))
t.summerNH<-t.summer.tmp[,which(lat>=0),]
t.winterSH<-t.summer.tmp[,which(lat<0),]
rm(nc)
t.winter<-abind(t.winterSH,t.winterNH,along=2)
t.summer<-abind(t.summerSH,t.summerNH,along=2)

nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/WOA/woa13_decav_t14_01v2.nc")
t.spring.tmp<-ncvar_get(nc,"t_an",start=c(1,1,1,1),count=c(-1,-1,40,-1))
t.springNH<-t.spring.tmp[,which(lat>=0),]
t.fallSH<-t.spring.tmp[,which(lat<0),]
rm(nc)
nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/WOA/woa13_decav_t16_01v2.nc")
t.fall.tmp<-ncvar_get(nc,"t_an",start=c(1,1,1,1),count=c(-1,-1,40,-1))
t.fallNH<-t.fall.tmp[,which(lat>=0),]
t.springSH<-t.fall.tmp[,which(lat<0),]
rm(nc)
t.spring<-abind(t.springSH,t.springNH,along=2)
t.fall<-abind(t.fallSH,t.fallNH,along=2)

#saturation state (NOTE: GLODAP has same latitude, but different longitude)
nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/GLODAP/GLODAPv2.2016b.OmegaC.nc")
OmegaC<-ncvar_get(nc,"OmegaC",start=c(1,1,1),count=c(-1,-1,-1))	
lon.glodap<-nc$dim$lon$vals
lon.glodap[which(lon.glodap>180)]<-lon.glodap[which(lon.glodap>180)]-360
depth.glodap<-ncvar_get(nc,"Depth",start=c(1),count=c(-1))	
rm(nc)
nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/GLODAP/GLODAPv2.2016b.pHtsinsitutp.nc")
pH<-ncvar_get(nc,"pHtsinsitutp",start=c(1,1,1),count=c(-1,-1,-1))	
rm(nc)

#***********************************************#

#coretop calibration of N.pachy.s
Npachy.d<-c(0,100)																									#assigned depth habitat
d<-which(depth<max(Npachy.d) & depth>min(Npachy.d))
d1<-which(depth.glodap<max(Npachy.d) & depth.glodap>min(Npachy.d))
Npachy.ta<-t.ann[,,d]
Npachy.ta<-apply(Npachy.ta,c(1,2),mean)
Npachy.t1<-t.summer[,,d]
Npachy.t1<-apply(Npachy.t1,c(1,2),mean)
Npachy.t2<-t.spring[,,d]
Npachy.t2<-apply(Npachy.t2,c(1,2),mean)
Npachy.s1<-s.summer[,,d]
Npachy.s1<-apply(Npachy.s1,c(1,2),mean)
Npachy.s2<-s.spring[,,d]
Npachy.s2<-apply(Npachy.s2,c(1,2),mean)
#Start with spring, but if mean annual T<=10, replace with summer based on Jonkers 15
Npachy.t<-Npachy.t2
Npachy.t[which(Npachy.ta<=10)]<-Npachy.t1[which(Npachy.ta<=10)]
Npachy.s<-Npachy.s2
Npachy.s[which(Npachy.ta<=10)]<-Npachy.s1[which(Npachy.ta<=10)]

Npachy.WOA.s<-seq(,length.out=nrow(MgCa.dat.Npachy))
Npachy.WOA.t<-seq(,length.out=nrow(MgCa.dat.Npachy))
Npachy.glodap.omega.s<-seq(,length.out=nrow(MgCa.dat.Npachy))
Npachy.glodap.omega<-seq(,length.out=nrow(MgCa.dat.Npachy))
Npachy.glodap.pH.s<-seq(,length.out=nrow(MgCa.dat.Npachy))
Npachy.glodap.pH<-seq(,length.out=nrow(MgCa.dat.Npachy))

for(i in 1:nrow(MgCa.dat.Npachy)) {
	lon.tmp<-which(abs(lon-MgCa.dat.Npachy$lon[i])==min(abs(lon-MgCa.dat.Npachy$lon[i])))[1]
	lat.tmp<-which(abs(lat-MgCa.dat.Npachy$lat[i])==min(abs(lat-MgCa.dat.Npachy$lat[i])))[1]
	if (is.na(Npachy.s[lon.tmp,lat.tmp])==TRUE) {																	#if there's no model data at the proxy site
			lon.tmp<-seq(lon.tmp-1,lon.tmp+1)																		#take the surrounding gridboxes
			lat.tmp<-seq(lat.tmp-1,lat.tmp+1)
			Npachy.WOA.s[i]<-mean(Npachy.s[lon.tmp,lat.tmp],na.rm=TRUE)
			Npachy.WOA.t[i]<-mean(Npachy.t[lon.tmp,lat.tmp],na.rm=TRUE)
		} else {
			Npachy.WOA.s[i]<-Npachy.s[lon.tmp,lat.tmp]
			Npachy.WOA.t[i]<-Npachy.t[lon.tmp,lat.tmp]
		}
	lon.tmp<-which(abs(lon.glodap-MgCa.dat.Npachy$lon[i])==min(abs(lon.glodap-MgCa.dat.Npachy$lon[i])))[1]
	lat.tmp<-which(abs(lat-MgCa.dat.Npachy$lat[i])==min(abs(lat-MgCa.dat.Npachy$lat[i])))[1]
	depth.tmp<-which(abs(depth.glodap-MgCa.dat.Npachy$depth[i])==min(abs(depth.glodap-MgCa.dat.Npachy$depth[i])))[1]
	if (is.na(OmegaC[lon.tmp,lat.tmp,depth.tmp])==TRUE) {																	#if there's no model data at the proxy site
			lon.tmp<-seq(lon.tmp-1,lon.tmp+1)																				#take the surrounding gridboxes
			lat.tmp<-seq(lat.tmp-1,lat.tmp+1)
			depth.tmp<-seq(depth.tmp-1,depth.tmp+1)
			Npachy.glodap.omega[i]<-mean(OmegaC[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE)
			Npachy.glodap.pH[i]<-mean(pH[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE)
			Npachy.glodap.omega.s[i]<-mean(OmegaC[lon.tmp,lat.tmp,d1],na.rm=TRUE)
			Npachy.glodap.pH.s[i]<-mean(pH[lon.tmp,lat.tmp,d1],na.rm=TRUE)
		} else {
			Npachy.glodap.omega[i]<-OmegaC[lon.tmp,lat.tmp,depth.tmp]
			Npachy.glodap.pH[i]<-pH[lon.tmp,lat.tmp,depth.tmp]
			Npachy.glodap.omega.s[i]<-mean(OmegaC[lon.tmp,lat.tmp,d1])
			Npachy.glodap.pH.s[i]<-mean(pH[lon.tmp,lat.tmp,d1])
		}
}

MgCa.dat.Npachy$WOA.s<-Npachy.WOA.s
MgCa.dat.Npachy$WOA.t<-Npachy.WOA.t
MgCa.dat.Npachy$GLODAP.OmegaC.surf<-Npachy.glodap.omega.s
MgCa.dat.Npachy$GLODAP.pH.surf<-Npachy.glodap.pH.s
MgCa.dat.Npachy$GLODAP.OmegaC.deep<-Npachy.glodap.omega
MgCa.dat.Npachy$GLODAP.pH.deep<-Npachy.glodap.pH

cut<-which(MgCa.dat.Npachy$GLODAP.OmegaC.deep<=OmegaCcut)
MgCa.dat.Npachy$GLODAP.OmegaC.deep.cut<-rep(OmegaCcut,length.out=nrow(MgCa.dat.Npachy))
MgCa.dat.Npachy$GLODAP.OmegaC.deep.cut[cut]<-MgCa.dat.Npachy$GLODAP.OmegaC.deep[cut]
MgCa.dat.Npachy$GLODAP.pH.deep.cut<-rep(mean(MgCa.dat.Npachy$GLODAP.pH.deep[-cut],na.rm=T),length.out=nrow(MgCa.dat.Npachy))
MgCa.dat.Npachy$GLODAP.pH.deep.cut[cut]<-MgCa.dat.Npachy$GLODAP.pH.deep[cut]

#***********************************************#
#coretop calibration of G.ruber
Gruber.d<-c(0,50)
d<-which(depth<max(Gruber.d) & depth>min(Gruber.d))
d1<-which(depth.glodap<max(Gruber.d) & depth.glodap>min(Gruber.d))
Gruber.t1<-t.ann[,,d]
Gruber.t1<-apply(Gruber.t1,c(1,2),mean)
Gruber.t2<-t.summer[,,d]
Gruber.t2<-apply(Gruber.t2,c(1,2),mean)
Gruber.s1<-s.ann[,,d]
Gruber.s1<-apply(Gruber.s1,c(1,2),mean)
Gruber.s2<-s.summer[,,d]
Gruber.s2<-apply(Gruber.s2,c(1,2),mean)
#if mean annual T<=25, replace with summer based on Jonkers 15
Gruber.t<-Gruber.t1
Gruber.t[which(Gruber.t1<=25)]<-Gruber.t2[which(Gruber.t1<=25)]
Gruber.s<-Gruber.s1
Gruber.s[which(Gruber.t1<=25)]<-Gruber.s2[which(Gruber.t1<=25)]

Gruber.WOA.s<-seq(,length.out=nrow(MgCa.dat.Gruber))
Gruber.WOA.t<-seq(,length.out=nrow(MgCa.dat.Gruber))
Gruber.glodap.omega.s<-seq(,length.out=nrow(MgCa.dat.Gruber))
Gruber.glodap.omega<-seq(,length.out=nrow(MgCa.dat.Gruber))
Gruber.glodap.pH.s<-seq(,length.out=nrow(MgCa.dat.Gruber))
Gruber.glodap.pH<-seq(,length.out=nrow(MgCa.dat.Gruber))

for(i in 1:nrow(MgCa.dat.Gruber)) {
	lon.tmp<-which(abs(lon-MgCa.dat.Gruber$lon[i])==min(abs(lon-MgCa.dat.Gruber$lon[i])))[1]
	lat.tmp<-which(abs(lat-MgCa.dat.Gruber$lat[i])==min(abs(lat-MgCa.dat.Gruber$lat[i])))[1]
	if (is.na(Gruber.s[lon.tmp,lat.tmp])==TRUE) {																	#if there's no model data at the proxy site
			lon.tmp<-seq(lon.tmp-1,lon.tmp+1)																		#take the surrounding gridboxes
			lat.tmp<-seq(lat.tmp-1,lat.tmp+1)
			Gruber.WOA.s[i]<-mean(Gruber.s[lon.tmp,lat.tmp],na.rm=TRUE)
			Gruber.WOA.t[i]<-mean(Gruber.t[lon.tmp,lat.tmp],na.rm=TRUE)
		} else {
			Gruber.WOA.s[i]<-Gruber.s[lon.tmp,lat.tmp]
			Gruber.WOA.t[i]<-Gruber.t[lon.tmp,lat.tmp]
		}
	lon.tmp<-which(abs(lon.glodap-MgCa.dat.Gruber$lon[i])==min(abs(lon.glodap-MgCa.dat.Gruber$lon[i])))[1]
	lat.tmp<-which(abs(lat-MgCa.dat.Gruber$lat[i])==min(abs(lat-MgCa.dat.Gruber$lat[i])))[1]
	depth.tmp<-which(abs(depth.glodap-MgCa.dat.Gruber$depth[i])==min(abs(depth.glodap-MgCa.dat.Gruber$depth[i])))[1]
	if (is.na(OmegaC[lon.tmp,lat.tmp,depth.tmp])==TRUE) {																	#if there's no model data at the proxy site
			lon.tmp<-seq(lon.tmp-1,lon.tmp+1)																				#take the surrounding gridboxes
			lat.tmp<-seq(lat.tmp-1,lat.tmp+1)
			depth.tmp<-seq(depth.tmp-1,depth.tmp+1)
			Gruber.glodap.omega[i]<-mean(OmegaC[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE)
			Gruber.glodap.pH[i]<-mean(pH[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE)
			Gruber.glodap.omega.s[i]<-mean(OmegaC[lon.tmp,lat.tmp,d1],na.rm=TRUE)
			Gruber.glodap.pH.s[i]<-mean(pH[lon.tmp,lat.tmp,d1],na.rm=TRUE)
		} else {
			Gruber.glodap.omega[i]<-OmegaC[lon.tmp,lat.tmp,depth.tmp]
			Gruber.glodap.pH[i]<-pH[lon.tmp,lat.tmp,depth.tmp]
			Gruber.glodap.omega.s[i]<-mean(OmegaC[lon.tmp,lat.tmp,d1])
			Gruber.glodap.pH.s[i]<-mean(pH[lon.tmp,lat.tmp,d1])
		}
}

MgCa.dat.Gruber$WOA.s<-Gruber.WOA.s
MgCa.dat.Gruber$WOA.t<-Gruber.WOA.t
MgCa.dat.Gruber$GLODAP.OmegaC.surf<-Gruber.glodap.omega.s
MgCa.dat.Gruber$GLODAP.pH.surf<-Gruber.glodap.pH.s
MgCa.dat.Gruber$GLODAP.OmegaC.deep<-Gruber.glodap.omega
MgCa.dat.Gruber$GLODAP.pH.deep<-Gruber.glodap.pH

cut<-which(MgCa.dat.Gruber$GLODAP.OmegaC.deep<=OmegaCcut)
MgCa.dat.Gruber$GLODAP.OmegaC.deep.cut<-rep(OmegaCcut,length.out=nrow(MgCa.dat.Gruber))
MgCa.dat.Gruber$GLODAP.OmegaC.deep.cut[cut]<-MgCa.dat.Gruber$GLODAP.OmegaC.deep[cut]
MgCa.dat.Gruber$GLODAP.pH.deep.cut<-rep(mean(MgCa.dat.Gruber$GLODAP.pH.deep[-cut]),length.out=nrow(MgCa.dat.Gruber))
MgCa.dat.Gruber$GLODAP.pH.deep.cut[cut]<-MgCa.dat.Gruber$GLODAP.pH.deep[cut]

#***********************************************#
#coretop calibration of G.inflata
Ginflata.d1<-c(74,151)											#75-150m north of 35; 250-350 south of 35º (including S. Atl). See Cleroux 07; Groenveld11.
Ginflata.d2<-c(249,351)
d1<-which(depth<max(Ginflata.d1) & depth>min(Ginflata.d1))
d2<-which(depth<max(Ginflata.d2) & depth>min(Ginflata.d2))
d11<-which(depth.glodap<max(Ginflata.d1) & depth.glodap>min(Ginflata.d1))
d21<-which(depth.glodap<max(Ginflata.d2) & depth.glodap>min(Ginflata.d2))

Ginflata.s1<-s.spring[,,d1]
Ginflata.s1<-apply(Ginflata.s1,c(1,2),mean)
Ginflata.t1<-t.spring[,,d1]
Ginflata.t1<-apply(Ginflata.t1,c(1,2),mean)
Ginflata.s2<-s.spring[,,d2]
Ginflata.s2<-apply(Ginflata.s2,c(1,2),mean)
Ginflata.t2<-t.spring[,,d2]
Ginflata.t2<-apply(Ginflata.t2,c(1,2),mean)

Ginflata.WOA.s<-seq(,length.out=nrow(MgCa.dat.Ginflata))
Ginflata.WOA.t<-seq(,length.out=nrow(MgCa.dat.Ginflata))
Ginflata.glodap.omega.s<-seq(,length.out=nrow(MgCa.dat.Ginflata))
Ginflata.glodap.omega<-seq(,length.out=nrow(MgCa.dat.Ginflata))
Ginflata.glodap.pH.s<-seq(,length.out=nrow(MgCa.dat.Ginflata))
Ginflata.glodap.pH<-seq(,length.out=nrow(MgCa.dat.Ginflata))

for(i in 1:nrow(MgCa.dat.Ginflata)) {
	lon.tmp<-which(abs(lon-MgCa.dat.Ginflata$lon[i])==min(abs(lon-MgCa.dat.Ginflata$lon[i])))[1]
	lat.tmp<-which(abs(lat-MgCa.dat.Ginflata$lat[i])==min(abs(lat-MgCa.dat.Ginflata$lat[i])))[1]
	if (lat[lat.tmp]>35) {
		Ginflata.s<-Ginflata.s1
		Ginflata.t<-Ginflata.t1
		} else {
		Ginflata.s<-Ginflata.s2
		Ginflata.t<-Ginflata.t2
		}
	if (is.na(Ginflata.s[lon.tmp,lat.tmp])==TRUE) {																	#if there's no model data at the proxy site
			lon.tmp<-seq(lon.tmp-1,lon.tmp+1)																		#take the surrounding gridboxes
			lat.tmp<-seq(lat.tmp-1,lat.tmp+1)
			Ginflata.WOA.s[i]<-mean(Ginflata.s[lon.tmp,lat.tmp],na.rm=TRUE)
			Ginflata.WOA.t[i]<-mean(Ginflata.t[lon.tmp,lat.tmp],na.rm=TRUE)
		} else {
			Ginflata.WOA.s[i]<-Ginflata.s[lon.tmp,lat.tmp]
			Ginflata.WOA.t[i]<-Ginflata.t[lon.tmp,lat.tmp]
		}
	lon.tmp<-which(abs(lon.glodap-MgCa.dat.Ginflata$lon[i])==min(abs(lon.glodap-MgCa.dat.Ginflata$lon[i])))[1]
	lat.tmp<-which(abs(lat-MgCa.dat.Ginflata$lat[i])==min(abs(lat-MgCa.dat.Ginflata$lat[i])))[1]
	depth.tmp<-which(abs(depth.glodap-MgCa.dat.Ginflata$depth[i])==min(abs(depth.glodap-MgCa.dat.Ginflata$depth[i])))[1]
	if (is.na(OmegaC[lon.tmp,lat.tmp,depth.tmp])==TRUE) {																	#if there's no model data at the proxy site
			lon.tmp<-seq(lon.tmp-1,lon.tmp+1)																				#take the surrounding gridboxes
			lat.tmp<-seq(lat.tmp-1,lat.tmp+1)
			depth.tmp<-seq(depth.tmp-1,depth.tmp+1)
			Ginflata.glodap.omega[i]<-mean(OmegaC[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE)
			Ginflata.glodap.pH[i]<-mean(pH[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE)
			if (lat[lat.tmp][2]>35) { 
				Ginflata.glodap.omega.s[i]<-mean(OmegaC[lon.tmp,lat.tmp,d11],na.rm=TRUE)
				Ginflata.glodap.pH.s[i]<-mean(pH[lon.tmp,lat.tmp,d11],na.rm=TRUE)
				} else {
				Ginflata.glodap.omega.s[i]<-mean(OmegaC[lon.tmp,lat.tmp,d21],na.rm=TRUE)
				Ginflata.glodap.pH.s[i]<-mean(pH[lon.tmp,lat.tmp,d21],na.rm=TRUE)	
				}
		} else {
			Ginflata.glodap.omega[i]<-OmegaC[lon.tmp,lat.tmp,depth.tmp]
			Ginflata.glodap.pH[i]<-pH[lon.tmp,lat.tmp,depth.tmp]
			if (lat[lat.tmp]>35) {
				Ginflata.glodap.omega.s[i]<-mean(OmegaC[lon.tmp,lat.tmp,d11],na.rm=TRUE)
				Ginflata.glodap.pH.s[i]<-mean(pH[lon.tmp,lat.tmp,d11],na.rm=TRUE)
				} else {
				Ginflata.glodap.omega.s[i]<-mean(OmegaC[lon.tmp,lat.tmp,d21],na.rm=TRUE)
				Ginflata.glodap.pH.s[i]<-mean(pH[lon.tmp,lat.tmp,d21],na.rm=TRUE)
				}
		}
}

MgCa.dat.Ginflata$WOA.s<-Ginflata.WOA.s
MgCa.dat.Ginflata$WOA.t<-Ginflata.WOA.t
MgCa.dat.Ginflata$GLODAP.OmegaC.surf<-Ginflata.glodap.omega.s
MgCa.dat.Ginflata$GLODAP.pH.surf<-Ginflata.glodap.pH.s
MgCa.dat.Ginflata$GLODAP.OmegaC.deep<-Ginflata.glodap.omega
MgCa.dat.Ginflata$GLODAP.pH.deep<-Ginflata.glodap.pH

cut<-which(MgCa.dat.Ginflata$GLODAP.OmegaC.deep<=OmegaCcut)
MgCa.dat.Ginflata$GLODAP.OmegaC.deep.cut<-rep(OmegaCcut,length.out=nrow(MgCa.dat.Ginflata))
MgCa.dat.Ginflata$GLODAP.OmegaC.deep.cut[cut]<-MgCa.dat.Ginflata$GLODAP.OmegaC.deep[cut]
MgCa.dat.Ginflata$GLODAP.pH.deep.cut<-rep(mean(MgCa.dat.Ginflata$GLODAP.pH.deep[-cut]),length.out=nrow(MgCa.dat.Ginflata))
MgCa.dat.Ginflata$GLODAP.pH.deep.cut[cut]<-MgCa.dat.Ginflata$GLODAP.pH.deep[cut]

#***********************************************#
#coretop calibration of G.bulloides
Gbulloides.d<-c(0,100)
d<-which(depth<max(Gbulloides.d) & depth>min(Gbulloides.d))
d1<-which(depth.glodap<max(Gbulloides.d) & depth.glodap>min(Gbulloides.d))
Gbulloides.ta<-t.ann[,,d]
Gbulloides.ta<-apply(Gbulloides.ta,c(1,2),mean)
Gbulloides.t1<-t.summer[,,d]
Gbulloides.t1<-apply(Gbulloides.t1,c(1,2),mean)
Gbulloides.t2<-t.spring[,,d]
Gbulloides.t2<-apply(Gbulloides.t2,c(1,2),mean)
Gbulloides.s1<-s.summer[,,d]
Gbulloides.s1<-apply(Gbulloides.s1,c(1,2),mean)
Gbulloides.s2<-s.spring[,,d]
Gbulloides.s2<-apply(Gbulloides.s2,c(1,2),mean)
#Start with spring, but if mean annual T<=10, replace with summer based on Jonkers 15
Gbulloides.t<-Gbulloides.t2
Gbulloides.t[which(Gbulloides.ta<=10)]<-Gbulloides.t1[which(Gbulloides.ta<=10)]
Gbulloides.s<-Gbulloides.s2
Gbulloides.s[which(Gbulloides.ta<=10)]<-Gbulloides.s1[which(Gbulloides.ta<=10)]

Gbulloides.WOA.s<-seq(,length.out=nrow(MgCa.dat.Gbulloides))
Gbulloides.WOA.t<-seq(,length.out=nrow(MgCa.dat.Gbulloides))
Gbulloides.glodap.omega.s<-seq(,length.out=nrow(MgCa.dat.Gbulloides))
Gbulloides.glodap.omega<-seq(,length.out=nrow(MgCa.dat.Gbulloides))
Gbulloides.glodap.pH.s<-seq(,length.out=nrow(MgCa.dat.Gbulloides))
Gbulloides.glodap.pH<-seq(,length.out=nrow(MgCa.dat.Gbulloides))

for(i in 1:nrow(MgCa.dat.Gbulloides)) {
	lon.tmp<-which(abs(lon-MgCa.dat.Gbulloides$lon[i])==min(abs(lon-MgCa.dat.Gbulloides$lon[i])))[1]
	lat.tmp<-which(abs(lat-MgCa.dat.Gbulloides$lat[i])==min(abs(lat-MgCa.dat.Gbulloides$lat[i])))[1]
	if (is.na(Gbulloides.s[lon.tmp,lat.tmp])==TRUE) {																	#if there's no model data at the proxy site
			lon.tmp<-seq(lon.tmp-1,lon.tmp+1)																		#take the surrounding gridboxes
			lat.tmp<-seq(lat.tmp-1,lat.tmp+1)
			Gbulloides.WOA.s[i]<-mean(Gbulloides.s[lon.tmp,lat.tmp],na.rm=TRUE)
			Gbulloides.WOA.t[i]<-mean(Gbulloides.t[lon.tmp,lat.tmp],na.rm=TRUE)
		} else {
			Gbulloides.WOA.s[i]<-Gbulloides.s[lon.tmp,lat.tmp]
			Gbulloides.WOA.t[i]<-Gbulloides.t[lon.tmp,lat.tmp]
		}
	lon.tmp<-which(abs(lon.glodap-MgCa.dat.Gbulloides$lon[i])==min(abs(lon.glodap-MgCa.dat.Gbulloides$lon[i])))[1]
	lat.tmp<-which(abs(lat-MgCa.dat.Gbulloides$lat[i])==min(abs(lat-MgCa.dat.Gbulloides$lat[i])))[1]
	depth.tmp<-which(abs(depth.glodap-MgCa.dat.Gbulloides$depth[i])==min(abs(depth.glodap-MgCa.dat.Gbulloides$depth[i])))[1]
	if (is.na(OmegaC[lon.tmp,lat.tmp,depth.tmp])==TRUE) {																	#if there's no model data at the proxy site
			lon.tmp<-seq(lon.tmp-1,lon.tmp+1)																				#take the surrounding gridboxes
			lat.tmp<-seq(lat.tmp-1,lat.tmp+1)
			depth.tmp<-seq(depth.tmp-1,depth.tmp+1)
			Gbulloides.glodap.omega[i]<-mean(OmegaC[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE)
			Gbulloides.glodap.pH[i]<-mean(pH[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE)
			Gbulloides.glodap.omega.s[i]<-mean(OmegaC[lon.tmp,lat.tmp,d1],na.rm=TRUE)
			Gbulloides.glodap.pH.s[i]<-mean(pH[lon.tmp,lat.tmp,d1],na.rm=TRUE)
		} else {
			Gbulloides.glodap.omega[i]<-OmegaC[lon.tmp,lat.tmp,depth.tmp]
			Gbulloides.glodap.pH[i]<-pH[lon.tmp,lat.tmp,depth.tmp]
			Gbulloides.glodap.omega.s[i]<-mean(OmegaC[lon.tmp,lat.tmp,d1])
			Gbulloides.glodap.pH.s[i]<-mean(pH[lon.tmp,lat.tmp,d1])
		}
}

MgCa.dat.Gbulloides$WOA.s<-Gbulloides.WOA.s
MgCa.dat.Gbulloides$WOA.t<-Gbulloides.WOA.t
MgCa.dat.Gbulloides$GLODAP.OmegaC.surf<-Gbulloides.glodap.omega.s
MgCa.dat.Gbulloides$GLODAP.pH.surf<-Gbulloides.glodap.pH.s
MgCa.dat.Gbulloides$GLODAP.OmegaC.deep<-Gbulloides.glodap.omega
MgCa.dat.Gbulloides$GLODAP.pH.deep<-Gbulloides.glodap.pH

cut<-which(MgCa.dat.Gbulloides$GLODAP.OmegaC.deep<=OmegaCcut)
MgCa.dat.Gbulloides$GLODAP.OmegaC.deep.cut<-rep(OmegaCcut,length.out=nrow(MgCa.dat.Gbulloides))
MgCa.dat.Gbulloides$GLODAP.OmegaC.deep.cut[cut]<-MgCa.dat.Gbulloides$GLODAP.OmegaC.deep[cut]
MgCa.dat.Gbulloides$GLODAP.pH.deep.cut<-rep(mean(MgCa.dat.Gbulloides$GLODAP.pH.deep[-cut]),length.out=nrow(MgCa.dat.Gbulloides))
MgCa.dat.Gbulloides$GLODAP.pH.deep.cut[cut]<-MgCa.dat.Gbulloides$GLODAP.pH.deep[cut]

#***********************************************#
#coretop calibration of all species
#all data, but isolate T dependence by accounting for confounding variables.
MgCa.dat.Npachy$Mg.Ca_adj<-MgCa.dat.Npachy$Mg.Ca-0.0016*MgCa.dat.Npachy$mean_size
MgCa.dat.Gruber$Mg.Ca_adj<-MgCa.dat.Gruber$Mg.Ca-2.7*MgCa.dat.Gruber$GLODAP.OmegaC.deep.cut+0.008*MgCa.dat.Gruber$mean_size
MgCa.dat.Ginflata$Mg.Ca_adj<-MgCa.dat.Ginflata$Mg.Ca-MgCa.dat.Ginflata$GLODAP.OmegaC.surf^-0.97
MgCa.dat.Gbulloides$Mg.Ca_adj<-MgCa.dat.Gbulloides$Mg.Ca*.58 #arbitrary fudge factor

MgCa.dat.all<-rbind(MgCa.dat.Npachy,MgCa.dat.Gruber,MgCa.dat.Ginflata,MgCa.dat.Gbulloides)

#N. pachyderma
f=function(temp,size,a,b1,b2) {a*exp(b1*temp)+b2*size}
tmp<-MgCa.dat.Npachy[which(is.na(MgCa.dat.Npachy$WOA.t)==FALSE & is.na(MgCa.dat.Npachy$mean_size)==FALSE),]
NpachyMultfit<-nls(Mg.Ca~f(WOA.t,mean_size,a,b1,b2),data=tmp,start=c(a=0.1,b1=0.1,b2=0.001))
NpachyMult_r2<-cor(tmp$Mg.Ca,predict(NpachyMultfit))^2
NpachyMult_RMSE<-sqrt(sum((predict(NpachyMultfit)-tmp$Mg.Ca)^2)/nrow(tmp))

#T-only
NpachyfitT<-lm(log(Mg.Ca)~WOA.t, data=tmp)
NpachyT_r2<-summary(NpachyfitT)$r.squared
NpachyT_RMSE<-sqrt(sum((exp(predict(NpachyfitT))-tmp$Mg.Ca)^2)/nrow(tmp))

#size-only
NpachyfitSZ<-lm(Mg.Ca~mean_size, data=tmp)
NpachyT_r2<-summary(NpachyfitSZ)$r.squared
NpachyT_RMSE<-sqrt(sum((exp(predict(NpachyfitSZ))-tmp$Mg.Ca)^2)/nrow(tmp))

#fixed T sensitivity 
f=function(temp,size,a1,b2) {a1*exp(0.081*temp)+b2*size}
NpachyMultfitFXD<-nls(Mg.Ca~f(WOA.t,mean_size,a1,b2),data=tmp,start=c(a1=0.4,b2=0.001))
NpachyMultFXD_r2<-cor(tmp$Mg.Ca,predict(NpachyMultfitFXD))^2
NpachyMultFXD_RMSE<-sqrt(sum((predict(NpachyMultfitFXD)-tmp$Mg.Ca)^2)/nrow(tmp))

dev.new(width=10, height=6)
par(mfrow=c(1,2))
plot(tmp$Mg.Ca,exp(predict(NpachyfitT)),pch=20,xlim=c(0,3),ylim=c(0,3),las=1,xlab="measured Mg/Ca (mmol/mol)",ylab="predicted Mg/Ca (mmol/mol)",main="T-only",cex=1.2)
abline(0,1)
M06<-which(tmp$Reference=="Meland et al. 2006")
Y08<-which(tmp$Reference=="Yu et al. 2008")
K09<-which(tmp$Reference=="Kozdon et al. 2009; Simstich et al. 2003")
V16<-which(tmp$Reference=="Vázquez Riveiros et al. 2016")
EG00<-which(tmp$Reference=="Elderfield and Ganssen 2000")
points(tmp$Mg.Ca[M06],exp(predict(NpachyfitT))[M06],pch=20,col="firebrick")
points(tmp$Mg.Ca[Y08],exp(predict(NpachyfitT))[Y08],pch=20,col="cornflowerblue")
points(tmp$Mg.Ca[K09],exp(predict(NpachyfitT))[K09],pch=20,col="gold3")
points(tmp$Mg.Ca[V16],exp(predict(NpachyfitT))[V16],pch=20,col="palegreen3")
points(tmp$Mg.Ca[EG00],exp(predict(NpachyfitT))[EG00],pch=20,col="deeppink4")

plot(tmp$Mg.Ca,predict(NpachyMultfit),pch=20,xlim=c(0,3),ylim=c(0,3),las=1,xlab="measured Mg/Ca (mmol/mol)",ylab="predicted Mg/Ca (mmol/mol)",main="T, size",cex=1.2)
abline(0,1)
points(tmp$Mg.Ca[M06],predict(NpachyMultfit)[M06],pch=20,col="firebrick")
points(tmp$Mg.Ca[Y08],predict(NpachyMultfit)[Y08],pch=20,col="cornflowerblue")
points(tmp$Mg.Ca[K09],predict(NpachyMultfit)[K09],pch=20,col="gold3")
points(tmp$Mg.Ca[V16],predict(NpachyMultfit)[V16],pch=20,col="palegreen3")
points(tmp$Mg.Ca[EG00],predict(NpachyMultfit)[EG00],pch=20,col="deeppink4")
legend("topleft",legend=c("M06","Y08","K09","V16","EG00","Other"), pch=rep(20,6), col=c("firebrick","cornflowerblue","gold3","palegreen3","deeppink4","black"),cex=0.8)

#G. ruber
f=function(temp,omega,size,a,b1,b2,b3) {a*exp(b1*temp)+b2*omega+b3*size}
tmp<-MgCa.dat.Gruber[which(is.na(MgCa.dat.Gruber$WOA.t)==FALSE & is.na(MgCa.dat.Gruber$mean_size)==FALSE & is.na(MgCa.dat.Gruber$GLODAP.OmegaC.deep)==FALSE),]
GruberMultfit<-nls(Mg.Ca~f(WOA.t,GLODAP.OmegaC.deep.cut, mean_size,a,b1,b2,b3),data=tmp,start=c(a=0.1,b1=0.1,b2=1,b3=0.001))
GruberMult_r2<-cor(tmp$Mg.Ca,predict(GruberMultfit))^2
GruberMult_RMSE<-sqrt(sum((predict(GruberMultfit)-tmp$Mg.Ca)^2)/nrow(tmp))

#T-only
GruberfitT<-lm(log(Mg.Ca)~WOA.t,data=tmp)
GruberT_RMSE<-sqrt(sum((exp(predict(GruberfitT))-tmp$Mg.Ca)^2)/nrow(tmp))

#Fixed T sensitivity
f=function(temp,omega,size,a1,b2,b3) {a1*exp(0.081*temp)+b2*omega+b3*size}
GruberMultfitFXD<-nls(Mg.Ca~f(WOA.t,GLODAP.OmegaC.deep.cut, mean_size,a1,b2,b3),data=tmp,start=c(a1=0.4,b2=1,b3=0.001))
GruberMultFXD_r2<-cor(tmp$Mg.Ca,predict(GruberMultfitFXD))^2
GruberMultFXD_RMSE<-sqrt(sum((predict(GruberMultfitFXD)-tmp$Mg.Ca)^2)/nrow(tmp))

dev.new(width=10, height=6)
par(mfrow=c(1,2))
plot(tmp$Mg.Ca,exp(predict(GruberfitT)),pch=20,xlim=c(1,7),ylim=c(1,7),las=1,xlab="measured Mg/Ca (mmol/mol)",ylab="predicted Mg/Ca (mmol/mol)",main="T-only",cex=1.2)
abline(0,1)
D02<-which(tmp$Reference=="Dekens et al. 2002")
MBB09<-which(tmp$Reference=="Mathien-Blard and Bassinot 2009")
A10<-which(tmp$Reference=="Arbuszewski et al. 2010")
K15<-which(tmp$Reference=="Khider et al. 2015")
C08<-which(tmp$Reference=="Cleroux et al. 2008")
S11<-which(tmp$Reference=="Sabbatini et al. 2011")
M11<-which(tmp$Reference=="Mohtadi et al. 2011")
X10<-which(tmp$Reference=="Xu et al., 2010")
H17<-which(tmp$Reference=="Hollstein et al., 2017")
EG00<-which(tmp$Reference=="Elderfield and Ganssen 2000")
points(tmp$Mg.Ca[D02],exp(predict(GruberfitT))[D02],pch=20,col="firebrick")
points(tmp$Mg.Ca[MBB09],exp(predict(GruberfitT))[MBB09],pch=20,col="cornflowerblue")
points(tmp$Mg.Ca[A10],exp(predict(GruberfitT))[A10],pch=20,col="gold3")
points(tmp$Mg.Ca[K15],exp(predict(GruberfitT))[K15],pch=20,col="palegreen3")
points(tmp$Mg.Ca[C08],exp(predict(GruberfitT))[C08],pch=20,col="cyan2")
points(tmp$Mg.Ca[S11],exp(predict(GruberfitT))[S11],pch=20,col="chocolate")
points(tmp$Mg.Ca[M11],exp(predict(GruberfitT))[M11],pch=20,col="grey50")
points(tmp$Mg.Ca[X10],exp(predict(GruberfitT))[X10],pch=20,col="darkseagreen4")
points(tmp$Mg.Ca[H17],exp(predict(GruberfitT))[H17],pch=20,col="blue4")
points(tmp$Mg.Ca[EG00],exp(predict(GruberfitT))[EG00],pch=20,col="deeppink4")
legend("topleft",legend=c("D02","MBB09","A10","K15","C08","S11","M11","X10","H17","EG00","Other"), pch=rep(20,11), col=c("firebrick","cornflowerblue","gold3","palegreen3","cyan2","chocolate","grey50","darkseagreen4","blue4","deeppink4","black"),cex=0.8)

plot(tmp$Mg.Ca,predict(GruberMultfit),pch=20,xlim=c(1,7),ylim=c(1,7),las=1,xlab="measured Mg/Ca (mmol/mol)",ylab="predicted Mg/Ca (mmol/mol)",main="T, Omega.deep,size",cex=1.2)
abline(0,1)
points(tmp$Mg.Ca[D02],predict(GruberMultfit)[D02],pch=20,col="firebrick")
points(tmp$Mg.Ca[MBB09],predict(GruberMultfit)[MBB09],pch=20,col="cornflowerblue")
points(tmp$Mg.Ca[A10],predict(GruberMultfit)[A10],pch=20,col="gold3")
points(tmp$Mg.Ca[K15],predict(GruberMultfit)[K15],pch=20,col="palegreen3")
points(tmp$Mg.Ca[C08],predict(GruberMultfit)[C08],pch=20,col="cyan2")
points(tmp$Mg.Ca[S11],predict(GruberMultfit)[S11],pch=20,col="chocolate")
points(tmp$Mg.Ca[M11],predict(GruberMultfit)[M11],pch=20,col="grey50")
points(tmp$Mg.Ca[X10],predict(GruberMultfit)[X10],pch=20,col="darkseagreen4")
points(tmp$Mg.Ca[H17],predict(GruberMultfit)[H17],pch=20,col="blue4")
points(tmp$Mg.Ca[EG00],predict(GruberMultfit)[EG00],pch=20,col="deeppink4")

#G. inflata
#T, pH.surf
f=function(temp,pH,a,b1,b2) {a*exp(b1*temp)+b2*pH}
tmp<-MgCa.dat.Ginflata[which(is.na(MgCa.dat.Ginflata$WOA.t)==FALSE & is.na(MgCa.dat.Ginflata$GLODAP.pH.surf)==FALSE),]
GinflataMultfit1<-nls(Mg.Ca~f(WOA.t,GLODAP.pH.surf,a,b1,b2),data=tmp,start=c(a=1,b1=0.1,b2=1))
GinflataMult1_r2<-cor(tmp$Mg.Ca,predict(GinflataMultfit1))^2
GinflataMult1_RMSE<-sqrt(sum((predict(GinflataMultfit1)-tmp$Mg.Ca)^2)/nrow(tmp))

#T, OmegaC.surf
GinflataMultfit2<-lm(log(Mg.Ca)~WOA.t+log(GLODAP.OmegaC.surf),data=tmp)
GinflataMultfit2_RMSE<-sqrt(sum((exp(predict(GinflataMultfit2))-tmp$Mg.Ca)^2)/nrow(tmp))

#T-only
GinflatafitT<-lm(log(Mg.Ca)~WOA.t,data=tmp)
GinflataT_RMSE<-sqrt(sum((exp(predict(GinflatafitT))-tmp$Mg.Ca)^2)/nrow(tmp))

#fixed T sensitivity 
f=function(temp,omega,a1,b2) {a1*exp(0.081*temp)*omega^b2}
GinflataMultfitFXD<-nls(Mg.Ca~f(WOA.t,GLODAP.OmegaC.surf,a1,b2),data=tmp,start=c(a1=0.4,b2=0.001))
GinflataMultFXD_r2<-cor(tmp$Mg.Ca,predict(GinflataMultfitFXD))^2
GinflataMultFXD_RMSE<-sqrt(sum((predict(GinflataMultfitFXD)-tmp$Mg.Ca)^2)/nrow(tmp))

dev.new(width=10, height=6)
par(mfrow=c(1,2))
plot(tmp$Mg.Ca,exp(predict(GinflatafitT)),pch=20,xlim=c(0,3),ylim=c(0,3),las=1,xlab="measured Mg/Ca (mmol/mol)",ylab="predicted Mg/Ca (mmol/mol)",main="T-only",cex=1.2)
abline(0,1)
C13<-which(tmp$Reference=="Cleroux et al. 2013")
C08<-which(tmp$Reference=="Cleroux et al. 2008")
Y08<-which(tmp$Reference=="Yu et al. 2008")
GC11<-which(tmp$Reference=="Groeneveld and Chiessi 2011")
EG00<-which(tmp$Reference=="Elderfield and Ganssen 2000")
points(tmp$Mg.Ca[C13],exp(predict(GinflatafitT))[C13],pch=20,col="firebrick")
points(tmp$Mg.Ca[Y08],exp(predict(GinflatafitT))[Y08],pch=20,col="cornflowerblue")
points(tmp$Mg.Ca[GC11],exp(predict(GinflatafitT))[GC11],pch=20,col="gold3")
points(tmp$Mg.Ca[C08],exp(predict(GinflatafitT))[C08],pch=20,col="cyan2")
points(tmp$Mg.Ca[EG00],exp(predict(GinflatafitT))[EG00],pch=20,col="deeppink4")
legend("topleft",legend=c("C13","C08","Y08","GC11","EG00","Other"), pch=rep(20,6), col=c("firebrick","cyan2","cornflowerblue","gold3","deeppink4","black"))

plot(tmp$Mg.Ca,exp(predict(GinflataMultfit2)),pch=20,xlim=c(0,3),ylim=c(0,3),las=1,xlab="measured Mg/Ca (mmol/mol)",ylab="predicted Mg/Ca (mmol/mol)",main="T, Omega.surf",cex=1.2)
abline(0,1)
points(tmp$Mg.Ca[C13],exp(predict(GinflataMultfit2))[C13],pch=20,col="firebrick")
points(tmp$Mg.Ca[Y08],exp(predict(GinflataMultfit2))[Y08],pch=20,col="cornflowerblue")
points(tmp$Mg.Ca[GC11],exp(predict(GinflataMultfit2))[GC11],pch=20,col="gold3")
points(tmp$Mg.Ca[C08],exp(predict(GinflataMultfit2))[C08],pch=20,col="cyan2")
points(tmp$Mg.Ca[EG00],exp(predict(GinflataMultfit2))[EG00],pch=20,col="deeppink4")

#G. bulloides
#T-only
GbulloidesfitT<-lm(log(Mg.Ca)~WOA.t,data=MgCa.dat.Gbulloides)
GbulloidesT_RMSE<-sqrt(sum((exp(predict(GbulloidesfitT))-MgCa.dat.Gbulloides$Mg.Ca)^2)/nrow(MgCa.dat.Gbulloides))

#Fixed T sensitivity
f=function(temp,a1) {a1*exp(0.081*temp)}
GbulloidesfitFXD<-nls(Mg.Ca~f(WOA.t,a1),data=MgCa.dat.Gbulloides,start=c(a1=0.4))
GbulloidesFXD_r2<-cor(MgCa.dat.Gbulloides$Mg.Ca,predict(GbulloidesfitFXD))^2
GbulloidesFXD_RMSE<-sqrt(sum((predict(GbulloidesfitFXD)-MgCa.dat.Gbulloides$Mg.Ca)^2)/nrow(MgCa.dat.Gbulloides))

dev.new(width=6, height=6)
plot(MgCa.dat.Gbulloides$Mg.Ca,exp(predict(GbulloidesfitT)),pch=20,xlim=c(0,9),ylim=c(0,9),las=1,xlab="measured Mg/Ca (mmol/mol)",ylab="predicted Mg/Ca (mmol/mol)",main="T-only",cex=1.2)
abline(0,1)
QK17<-which(MgCa.dat.Gbulloides$Reference=="Quintana Krupinski 2017")
Y08<-which(MgCa.dat.Gbulloides$Reference=="Yu et al. 2008")
C08<-which(MgCa.dat.Gbulloides$Reference=="Cleroux et al. 2008")
V16<-which(MgCa.dat.Gbulloides$Reference=="Vázquez Riveiros et al. 2016")
M11<-which(MgCa.dat.Gbulloides$Reference=="Mohtadi et al. 2011")
EG00<-which(MgCa.dat.Gbulloides$Reference=="Elderfield and Ganssen 2000")
points(MgCa.dat.Gbulloides$Mg.Ca[QK17],exp(predict(GbulloidesfitT))[QK17],pch=20,col="gold3")
points(MgCa.dat.Gbulloides$Mg.Ca[Y08],exp(predict(GbulloidesfitT))[Y08],pch=20,col="cornflowerblue")
points(MgCa.dat.Gbulloides$Mg.Ca[C08],exp(predict(GbulloidesfitT))[C08],pch=20,col="cyan2")
points(MgCa.dat.Gbulloides$Mg.Ca[V16],exp(predict(GbulloidesfitT))[V16],pch=20,col="palegreen3")
points(MgCa.dat.Gbulloides$Mg.Ca[M11],exp(predict(GbulloidesfitT))[M11],pch=20,col="grey50")
points(MgCa.dat.Gbulloides$Mg.Ca[EG00],exp(predict(GbulloidesfitT))[EG00],pch=20,col="deeppink4")
legend("topleft",legend=c("QK17","C08","Y08","V16","M11","EG00","Other"), pch=rep(20,6), col=c("gold3","cyan2","cornflowerblue","palegreen3","grey50","deeppink4","black"))

#All species
#T-only
tmp<-MgCa.dat.all[which(is.na(MgCa.dat.all$WOA.t)==FALSE),]
allfitT<-lm(log(Mg.Ca)~WOA.t,data=tmp)
allT_RMSE<-sqrt(sum((exp(predict(allfitT))-tmp$Mg.Ca)^2)/nrow(tmp))

#T-only with adjusted MgCa (without G bulloides) 
tmp<-MgCa.dat.all[which(is.na(MgCa.dat.all$WOA.t)==FALSE & is.na(MgCa.dat.all$Mg.Ca_adj)==FALSE),]
allfitT_adj<-lm(log(Mg.Ca_adj)[1:724]~WOA.t[1:724],data=tmp)
allT_adj_RMSE<-sqrt(sum((exp(predict(allfitT_adj))[1:724]-tmp$Mg.Ca_adj[1:724])^2)/724)

#T-only with adjusted MgCa (all) 
allfitT_adj2<-lm(log(Mg.Ca_adj)~WOA.t,data=tmp)
allT_adj2_RMSE<-sqrt(sum((exp(predict(allfitT_adj2))-tmp$Mg.Ca_adj)^2)/nrow(tmp))

dev.new(width=8, height=6)
par(mfrow=c(1,2))
plot(MgCa.dat.all$WOA.t[1:724],MgCa.dat.all$Mg.Ca_adj[1:724],pch=20,xlim=c(-2,32),ylim=c(0,9),las=1,xlab="T (ºC)",ylab="Mg/Ca (mmol/mol)",main="minus secondary",cex=1.2)
legend("topleft",legend=c("N. pachyderma","G. ruber","G. inflata","G. bulloides"), pch=rep(20,4), col=c("palegreen3","cornflowerblue","firebrick","gold3"),bty="n")
points(MgCa.dat.Npachy$WOA.t,MgCa.dat.Npachy$Mg.Ca_adj,pch=20,col="palegreen3")
points(MgCa.dat.Gruber$WOA.t,MgCa.dat.Gruber$Mg.Ca_adj,pch=20,col="cornflowerblue")
points(MgCa.dat.Ginflata$WOA.t,MgCa.dat.Ginflata$Mg.Ca_adj,pch=20,col="firebrick")
points(MgCa.dat.Gbulloides$WOA.t,MgCa.dat.Gbulloides$Mg.Ca,pch=20,cex=1.2)
points(MgCa.dat.Gbulloides$WOA.t,MgCa.dat.Gbulloides$Mg.Ca,pch=20,col="gold3")
lines(seq(-4,30),exp(allfitT_adj$coeff[1])*exp(seq(-4,30)*allfitT_adj$coeff[2]))

plot(MgCa.dat.all$WOA.t,MgCa.dat.all$Mg.Ca_adj,pch=20,xlim=c(-2,32),ylim=c(0,9),las=1,xlab="T (ºC)",ylab="Mg/Ca (mmol/mol)",main="empirical G.bul",cex=1.2)
points(MgCa.dat.Npachy$WOA.t,MgCa.dat.Npachy$Mg.Ca_adj,pch=20,col="palegreen3")
points(MgCa.dat.Gruber$WOA.t,MgCa.dat.Gruber$Mg.Ca_adj,pch=20,col="cornflowerblue")
points(MgCa.dat.Ginflata$WOA.t,MgCa.dat.Ginflata$Mg.Ca_adj,pch=20,col="firebrick")
points(MgCa.dat.Gbulloides$WOA.t,MgCa.dat.Gbulloides$Mg.Ca_adj,pch=20,col="gold3")
lines(seq(-4,30),exp(allfitT_adj$coeff[1])*exp(seq(-4,30)*allfitT_adj$coeff[2]))