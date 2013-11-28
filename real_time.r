

library(fields)
library(rworldmap)
library(sp)
library(raster)
library(maptools)
library(gstat)
library(car)
library(stringr)


## All Oz Weather

# Read the data as fixed-width format
in_data=read.fwf("ftp://ftp2.bom.gov.au/anon/gen/fwo/IDY03021.txt",
                 widths=c(12,5,6,3,5,4,2,3,1,4,2,4,1,7,3,1,3,7,3,2),
                 col.names=c("Station",
                             "Lat",
                             "Lon",
                             "Day",
                             "Hour",
                             "Vis",
                             "Cld",
                             "WDir",
                             "S1",
                             "WSpd",
                             "Temp",
                             "RH",
                             "BQ",
                             "Pres",
                             "Rain",
                             "S2",
                             "RHour",
                             "Weather",
                             "Tx",
                             "Tn"),
                 skip=2,stringsAsFactors=F)

# Remove all the state data file identifiers
in_data=subset(in_data,substr(in_data$Station,1,3)!="IDY")
# Remove all the metadata lines
in_data=subset(in_data,substr(in_data$Station,1,1)!="|")
in_data=subset(in_data,substr(in_data$Station,1,1)!="+")

# We don't care about sites without observations
in_data=subset(in_data,in_data$Temp != "-- ")
in_data=subset(in_data,in_data$RH != "-- ")
in_data=subset(in_data,in_data$WSpd != "--- ")


# We don't care about Antarctica or anything after
ant=which(in_data$Station=="Casey       ")
in_data=in_data[1:(ant-1),]

# Some stations marked as "Calm" - give them a wind speed instead
in_data$WSpd[in_data$WDir==" Ca"]=0
in_data$WDir[in_data$WDir==" Ca"]="Calm"

WDirVec=c("N","NNE","NE","ENE","E","ESE","SE","SSE","S","SSW","SW","WSW","W","WNW","NW","NNW","Calm")
WDirAng=c(seq(0,359,by=22.5),0)
rcstring=paste("'",WDirVec,"'=",WDirAng,sep="",collapse=";")
in_data$WindAng=recode(str_trim(in_data$WDir),rcstring)





# Make the fields we're interested in numeric (and metric)
in_data$Temp = as.numeric(in_data$Temp)
in_data$RH = as.numeric(in_data$RH)
in_data$WSpd = as.numeric(in_data$WSpd)*1.852


# Convert latitudes and longitudes to decimal
in_data$Lat=-(as.numeric(substr(in_data$Lat,1,2))+(0.01*as.numeric((substr(in_data$Lat,3,4)))))
in_data$Lon=(as.numeric(substr(in_data$Lon,1,3))+(0.01*as.numeric((substr(in_data$Lon,4,5)))))


# Ignore eastern islands
in_data=subset(in_data,in_data$Lon < 155)

# Data set contains old data, from stations with infrequent reporting.  We only want reports
# From the last few hours
in_data$aHour = substr(in_data$Hour,1,2)
hourF=table(in_data$aHour)
hourS=rev(sort(hourF))[1:2]
sel_hour=names(hourS)
in_data=subset(in_data,in_data$aHour %in% sel_hour)

in_data=in_data[complete.cases(in_data),]

# Duplicate stations are present - something to do with rainfall pre/post 9am?
unq=!duplicated(in_data$Station)
in_data=in_data[unq,]


DF=raster("DF.tif")
FFDI=function(Temp,RH,Wind,DF){
  F <- 2 * exp(-.45 + .987 * log(DF + .001) - .0345 * RH + .0338 * Temp + .0234 * Wind) 
  F
}

ele=raster("DEM.tif")


in_data$x = in_data$Lon
in_data$y = in_data$Lat
coordinates(in_data)=c("Lon","Lat")
in_data$ele = extract(ele,in_data)
in_data$ele[is.na(in_data$ele)]=5

proj4string(in_data)="+proj=longlat +datum=WGS84"
in_data=remove.duplicates(in_data)

xy=data.frame(xyFromCell(ele,1:ncell(ele)))
xy$ele <- values(ele)
v=variogram(Temp ~ ele,in_data)
m<-fit.variogram(v,vgm(1,"Sph",200,1))
kr=gstat(NULL,"Temp",Temp~ele,data=in_data,model=m)
rTemp=ele
coordinates(xy)=c("x","y")
proj4string(xy)="+proj=longlat +datum=WGS84"
kr.p=predict(kr,newdata=xy)
proj4string(rTemp)="+proj=longlat +datum=WGS84"
values(rTemp)=as.vector(kr.p@data$Temp.pred)

xy=data.frame(xyFromCell(ele,1:ncell(ele)))
xy$ele <- values(ele)
v=variogram(WSpd ~ ele,in_data)
m<-fit.variogram(v,vgm(1,"Sph",200,1))
kr=gstat(NULL,"WSpd",WSpd~ele,data=in_data,model=m)
rWind=ele
coordinates(xy)=c("x","y")
proj4string(xy)="+proj=longlat +datum=WGS84"
kr.p=predict(kr,newdata=xy)
proj4string(rWind)="+proj=longlat +datum=WGS84"
values(rWind)=as.vector(kr.p@data$WSpd.pred)

xy=data.frame(xyFromCell(ele,1:ncell(ele)))
xy$ele <- values(ele)
v=variogram(RH ~ ele,in_data)
m<-fit.variogram(v,vgm(1,"Sph",200,1))
kr=gstat(NULL,"RH",RH~ele,data=in_data,model=m)
rRH=ele
coordinates(xy)=c("x","y")
proj4string(xy)="+proj=longlat +datum=WGS84"
kr.p=predict(kr,newdata=xy)
proj4string(rRH)="+proj=longlat +datum=WGS84"
values(rRH)=as.vector(kr.p@data$RH.pred)

rDF = ele
proj4string(rDF)="+proj=longlat +datum=WGS84"
values(rDF)=7

rFFDI=FFDI(rTemp,rRH,rWind,rDF)


#### Define official FFDI category Palette;
pal_low=colorRampPalette(c("#5d9450","#74c783"))(12)
pal_high=colorRampPalette(c("#317c78","#3697c2"))(13)
pal_vhigh=colorRampPalette(c("#9cb057","#f3d948"))(25)
pal_severe=colorRampPalette(c("#a37515","#e88333"))(25)
pal_extreme=colorRampPalette(c("#890e00","#e40c01"))(25)
pal_cat=colorRampPalette(c("#514657","#7d7d7d"))(100)
all_pal=c(pal_low,pal_high,pal_vhigh,pal_severe,pal_extreme,pal_cat)

aus=readShapePoly("aus.shp")
png("ffdi_map.png",width=1000,height=1000)
par(mai=c(1,1,1,1))
plot(rFFDI,col=all_pal,main=paste("FFDI: ",Sys.time(),sep=""),zlim=c(0,200))
plot(aus,add=T,lwd=2)
dev.off()

