library(sp)
library(maptools)
library(lubridate)
library (openxlsx)
library (gtools)
library (geosphere)
library(reshape2)
library (tidyverse)
library(readxl)
library (gifski)
library(viridis)
library (sf)
library(GISTools)
library(rgdal)
library (gganimate)
library(raster)
library(rgeos)
library(dplyr)
dat <- read.csv(".\\data\\nw-wrs-data.csv", 
                stringsAsFactors=FALSE, header=T, fileEncoding="latin1")
load("data\\stratum_map.RData")
load("data\\bcr_map.Rdata")
load("data\\lcc_map.Rdata")
load("data\\state_prov_map.Rdata")
load("data\\circle_2013_map.RData")
stratmod <- readOGR("C:\\Users\\rolek.brian\\Documents\\Projects\\northwest-wrs\\data\\shapefiles\\stratum_map_modified.shp")

########################
# join spatial information 
# to designate strata for each 
# survey
########################
# add up BAEG columns
dat$BAEG_all <- dat$BAEG_noage + dat$BAEGad + dat$BAEGsa

# Fix dates
for (i in 1:nrow(dat)){
  a <- gregexpr("-", dat$Date[i])[[1]]
  b <- substr(dat$Date[i],a[[1]][1]+1, nchar(dat$Date[i]))
  if (!is.na(as.numeric(dat$Date[i]))) {
    if (as.numeric(dat$Date[i])>30000){ #[950]
    dat$Date2[i] <- format(convertToDate(dat$Date[i]), "%m-%d-%y")
    next
    }} else{
    if (is.na(dat$Date[i])) { dat$Date2[i] <- NA; next}
    if ( length(a)==1 & a[[1]]!=-1){
      dat$Date2[i] <- paste( substr(dat$Date[i],1,a[[1]][1]-1),
                           "-01-", ifelse(b>50,"19","20"),
                           substr(dat$Date[i],a[[1]][1]+1, nchar(dat$Date[i])),
                           sep="")
  }
    if ( length(a)==2){ dat$Date2[i] <- dat$Date[i] }
    if ( (!length(a) %in% c(1,2)) | 
       a[[1]]==-1 | 
       any(gregexpr("&", dat$Date[i])[[1]]!=-1) |
       any(gregexpr(",", dat$Date[i])[[1]]!=-1) )  
  { dat$Date2[i] <- NA }
}}

dat$Date_form <- NA
for (i in 1:dim(dat)[[1]]){
  if (is.na(dat$Date2[i])) { next } else{
    len <- nchar(dat$Date2[i])-gregexpr("-", dat$Date2[i])[[1]][2]
    if (len==4){
      dat$Date_form[i] <- as.character(as.POSIXlt(dat$Date2[i], format="%m-%d-%Y"))
    } else{
      if (len==2){
        dat$Date_form[i] <- as.character(as.POSIXlt(dat$Date2[i], format="%m-%d-%y"))
      } else{ next }
    }}}

# manipulate dates
dat$Date_form <- as.POSIXlt(dat$Date_form, format=c("%Y-%m-%d") )
dat$month <- month(dat$Date_form)
dat$year <- year(dat$Date_form)
dat$day <- dat$Date_form$yday
write.csv(dat, file="C:/Users/rolek.brian/Documents/Projects/northwest-wrs/data/nw-wrs-data_forLeah.csv")

dat$year_month <- paste(dat$year, "_", dat$month, sep="")
dat <- dat[order(dat$route, dat$Date_form),]
dat <- dat[dat$year<=2020,]

dat$lat <- as.numeric(dat$lat)
dat$long <- as.numeric(dat$long)
dat <- dat[dat$long>= -180,]
dat <- dat[ ! (is.na(dat$lat) | is.na(dat$long)), ]
dat <- dat[, c(1:8, 50, 9:56)]
##
# subset to dates between 
# 14 Dec and 05 Jan for consistency 
# with CBC data
dat$dayOfMonth <- day(dat$Date_form) 
# remove nonsurveyed sites where 
# day of month was not provided
dat <- dat[nchar(dat$Date)>6,]
# Subset by date, 14 Dec and 05 Jan
dat <- dat[  ( (dat$month==12 & dat$dayOfMonth>=14) & 
                 (dat$month==12 & dat$dayOfMonth<=31) ) |
               ( (dat$month==1 & dat$dayOfMonth<=5)  &
                   (dat$month==1 & dat$dayOfMonth>=1) ), ]

uniq <- dat[!duplicated(dat[, c(49,50)]),]
spdf <- SpatialPointsDataFrame(coords = data.frame(uniq$long, uniq$lat ), data = uniq,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
spdf<- spTransform(spdf,crs(stratmod))
strat<- over(spdf, stratmod, returnList = F)
all <- cbind(spdf@data, strat, coordinates(spdf))
# for each survey
spdf.dat <- SpatialPointsDataFrame(coords = data.frame(dat$long, dat$lat ), data = dat,
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
spdf.dat<- spTransform(spdf.dat,crs(stratmod))
strat.dat<- over(spdf.dat, stratmod, returnList = F)
all.dat <- cbind(spdf.dat@data, strat.dat, coordinates(spdf.dat))
# for CBC surveys
circles<- spTransform(circle_2013_map, crs(stratmod))
cbc_strat <- over(circles, stratmod, returnList = F)
cbc <- cbind(circles@data, cbc_strat, coordinates(circles))


###################
# Summary stats
###################
# number of detections, mean num of detections per survey
samps <- merge(
  colSums(all.dat[, c(6:45)], na.rm=T),
  round(colMeans(all.dat[, c(6:45)], na.rm=T),3),
  by=0
)
samps[order(samps[,2], decreasing=T),]
# number of surveys total
all.dat.surveyed <- all.dat[!is.na(all.dat$total),]
dim(all.dat.surveyed)[[1]]
# number of routes surveyed
length(unique(all.dat.surveyed$route))
tab.survs.route.yrmo <- tapply(all.dat.surveyed$route, list(all.dat.surveyed$route, all.dat.surveyed$year_month) , length, default=NA)
tab.survs.route <- tapply(all.dat.surveyed$route, list(all.dat.surveyed$route, all.dat.surveyed$year) , length, default=0)
tab.survs.route
# num years a route was surveyed
tot.survs.route <- rowSums(tab.survs.route>0)
sum(tot.survs.route<=1)
sum(tot.survs.route<=7)
# routes surveyed per strat and year
tapply(all.dat.surveyed$unit_code , list(all.dat.surveyed$unit_code, all.dat.surveyed$year) , length, default=0)
tab.survs.unit <- tapply(all.dat.surveyed$unit_code , list(all.dat.surveyed$unit_code, all.dat.surveyed$year_month) , length, default=0)
tab.survs.unit <- tab.survs.unit[, mixedsort(colnames(tab.survs.unit))]
tab.survs.unit
# omit non surveyed data entries
all.dat <- all.dat[!is.na(all.dat$total),]
# add NA surveys so all are included
mtab <- melt(tab.survs.route.yrmo)
colnames(mtab) <- c("route", "year_month", "count")
mtab <- mtab[is.na(mtab$count),]
df.filler <- data.frame(X=NA, Date=NA, Miles=mean(all.dat$Miles, na.rm=T), Hours=NA, 
                        total=NA, RTHA=NA, AMKE=NA, NOHA=NA, 
                        BAEG_all=NA, BAEGad=NA, BAEGsa=NA, BAEG_noage=NA,
                        GOEA=NA, EAGLE=NA, RLHA=NA, RSHA=NA, FEHA=NA, SWHA=NA, BUTEO=NA,
                        WTKI=NA, PEFA=NA, PRFA=NA, MERL=NA, GYRF=NA, FALCON=NA,
                        COHA=NA, SSHA=NA, NOGO=NA, ACCIP=NA, OSPR=NA, GHOW=NA,
                        BNOW=NA, BUOW=NA, SEOW=NA, WESO=NA, NOPO=NA, LEOW=NA,
                        NOSO=NA, NOHO=NA, SNOW=NA, GGOW=NA, RAPTOR=NA, BADO=NA,
                        OWL=NA, ZTHA=NA, route=mtab$route,  
                        BAEG_all.1=NA, Date2=NA, Date_form=NA,
                        month=NA, year=NA, day=NA, year_month=mtab$year_month, dayOfMonth=NA)
filler.merge <- merge(spdf@data[,c(46,47,48,49,50)], df.filler, by="route", all.y=T)
filler.merge$lat[filler.merge$route=="BURLEY SW (ID)"] <- all.dat[all.dat$route=="BURLEY SW (ID)",]$lat[1]
filler.merge$long[filler.merge$route=="BURLEY SW (ID)"] <- all.dat[all.dat$route=="BURLEY SW (ID)",]$long[1]
filler.merge$lat[filler.merge$route=="PUGET ISLAND (CATHLAMET) (WA)"] <- all.dat[all.dat$route=="PUGET ISLAND (CATHLAMET) (WA)",]$lat[1]
filler.merge$long[filler.merge$route=="PUGET ISLAND (CATHLAMET) (WA)"] <- all.dat[all.dat$route=="PUGET ISLAND (CATHLAMET) (WA)",]$long[1]

spdf.filler <- SpatialPointsDataFrame(coords = data.frame(filler.merge$long, filler.merge$lat ), data = filler.merge,
                                      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
spdf.filler <- spTransform(spdf.filler,crs(stratmod))
strat.filler <- over(spdf.filler, stratmod, returnList = F)
all.filler <- cbind(spdf.filler@data, strat.filler, coordinates(spdf.filler))
all.filler$year <- as.numeric(substr(all.filler$year_month,1,4))
for (i in 1:dim(all.filler)[[1]]){
  a <- gregexpr("_", all.filler$year_month[i])[[1]]
  all.filler$month[i] <- as.numeric(substr(all.filler$year_month[i],a[[1]][1]+1, nchar(as.character(all.filler$year_month[i])))) 
}
all.filler$Date <- all.filler$Date2 <- 
  paste(all.filler$year, "-", all.filler$month, "-01", sep="")
all.filler$Date_form <- as.POSIXlt(all.filler$Date2, format="%Y-%m-%d")
all.filler$day <- all.filler$Date_form$yday
all.filler <- all.filler[, c(6:50,2,1,3:5,51:65)]
colnames(all.filler)[c(64,65)] <- c("dat.long", "dat.lat")
all.dat <- rbind (all.filler, all.dat)
# remove California 
all.dat <- all.dat[all.dat$state!="CA",] 
all.dat <- all.dat[!is.na(all.dat$route),]
all.dat$state_num <- as.numeric(factor(all.dat$state))
# remove routes where centroids fall outside of strata
all.dat <- all.dat[!is.na(all.dat$bcr_num),]
# subset data to routes with >1 survey
# num years a route was surveyed
keeper.routes <- rownames(tab.survs.route>1)
all.dat <- all.dat[all.dat$route %in% keeper.routes,]
all.dat <- all.dat[!is.na(all.dat$route),]
all.dat <- all.dat[!all.dat$month %in% c(11,4),  ]
all.dat <- all.dat[!is.na(all.dat$Miles), ]
ym <- unique(all.dat$year_month)
surv_yr<- c("04-05", "04-05", "05-06", "05-06",
            "06-07", "06-07", "07-08", "07-08",
            "08-09", "08-09", "09-10", "09-10",
            "10-11", "10-11", "11-12", "11-12",
            "12-13", "12-13", "13-14", "13-14",
            "14-15", "14-15", "15-16", "15-16",
            "16-17", "16-17", "17-18", "17-18",
            "18-19", "18-19", "19-20", "19-20")
uniq.time <- data.frame( ym=mixedsort(ym),
                         y=as.numeric(substr(ym,1,4)),
                         m=NA, surv_yr=surv_yr)
all.dat <- merge(all.dat, uniq.time[,c("ym", "surv_yr")], 
                 by.x="year_month", by.y="ym")
all.dat$surv_yr <- factor(all.dat$surv_yr, levels = unique(uniq.time$surv_yr) )
# Calculate speed km/h, impute average for values that are too fast >100 kph
all.dat$sp_kph <- (all.dat$Miles*1.60934)/all.dat$Hours
all.dat$sp_std <- (all.dat$sp_kph-mean(all.dat$sp_kph, na.rm=T))/sd(all.dat$sp_kph, na.rm=T)
all.dat$sp_std[all.dat$sp_kph>=100 | is.na(all.dat$sp_kph)] <- 0
##
detections <- sort(colSums(all.dat[, 7:46], na.rm=T), decreasing=T)
keeper.spp <- detections[detections>=500]
keeper.spp <- names(keeper.spp[!names(keeper.spp) %in% c("BAEGad", "BAEGsa", "BUTEO", "BAEG_noage")])
keeper.spp<- keeper.spp[-16]

all.dat$time_num <- as.numeric(all.dat$surv_yr)
all.dat$route_num <- as.numeric(factor(all.dat$route))

# make an index to average routes over strata
all.dat$strata_num <- as.numeric(factor(all.dat$unit_code))
all.dat$bcr_num2 <- as.numeric(factor(all.dat$bcr_num))
all.dat$lcc_num2 <- as.numeric(factor(all.dat$lcc_num))
all.dat$state_num <- as.numeric(factor(all.dat$unit_name))
surv.dat <- all.dat[!is.na(all.dat$total),]
ind <- table(surv.dat$strata_num, surv.dat$time_num, surv.dat$route_num)

# fill in NAs for BCR and LCC
routes  <- all.dat[!duplicated(all.dat$route),]
routes <- routes[order(routes$route_num), ]
ind.bcr <- table(routes$bcr_num2)
ind.lcc <- table(routes$lcc_num2)
ind.strata <- table(routes$strata_num)
ind.state <- table(routes$state_num)
bcr <- array(NA, dim=c(length(ind.bcr), max(ind.bcr)) )
lcc <- array(NA, dim=c(length(ind.lcc), max(ind.lcc)) )
strata <- array(NA, dim=c(length(ind.strata), max(ind.strata)) )
state <- array(NA, dim=c(length(ind.state), max(ind.state)) )

for (z in 1:3){
  w <- which(routes$bcr_num2==z)
  bcr[z,1:length(w)] <- w
  w2 <- which(routes$lcc_num2==z)
  lcc[z,1:length(w2)] <- w2
}
for (zz in 1:length(ind.strata)){
  w <- which(routes$strata_num==zz)
  strata[zz,1:length(w)] <- w
}
for (zzz in 1:length(ind.state)){
  w <- which(routes$state_num==zzz)
  state[zzz,1:length(w)] <- w
}

#################
# set up data for JAGS
###################
# remove NA routes
all.dat <- all.dat[!is.na(all.dat$route),]
# save spatial data file
wrs.sp <- SpatialPointsDataFrame(coords = data.frame(all.dat$long, all.dat$lat ), data = all.dat,
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
wrs.sp<- spTransform(wrs.sp,crs(stratmod))
writeOGR(obj=wrs.sp, dsn="C://Users//rolek.brian//Documents//Projects//northwest-wrs//data//shapefiles", 
         layer="wrs counts", driver="ESRI Shapefile")

ind2 <- table( surv.dat$route_num, surv.dat$time_num)
which.fun <- function(x){  min(which(x>0)) }
first.surv <- apply(ind2, 1, which.fun)
all.dat <- merge(all.dat, first.surv, by.x="route_num", by.y=0 )
all.dat$first.surv.yr <- all.dat$y+2003
all.dat99 <- all.dat[!is.na(all.dat$total),]
frst.wrs <- all.dat99[(all.dat99$time_num+2003)==all.dat99$first.surv.yr,]
frst.wrs <- frst.wrs[!duplicated(frst.wrs$route_num), ]
dwrs <- list()
spp.list <- temp <-  c()
# calculate first survey for each route

all.dat99<- all.dat99[order(all.dat99$route_num, all.dat99$time_num), ]
#frst.wrs <- all.dat99[!duplicated(all.dat99$route_num), ]$time_num
#names(frst.wrs) <- all.dat99[!duplicated(all.dat99$route_num), ]$route_num

for (sp in 1:length(keeper.spp)){
dwrs[[sp]] <- list( 
    count=all.dat[ ,keeper.spp[sp] ],
    dist=(all.dat$Miles*1.609)/100, # per 100 km
    time=all.dat$time_num,
    rt=all.dat$route_num,
    frst= frst.wrs$y ,
    mn1st= frst.wrs[, keeper.spp[sp]],
    str_ind= strata,
    nroutesInStrata= ind.strata, 
    state_ind=state,
    nroutesInState=ind.state,
    bcr_ind=bcr,
    nroutesInBCR=ind.bcr,
    lcc_ind=lcc,
    nroutesInLCC=ind.lcc, 
    sp=all.dat$sp_std, 
    ncounts=dim(all.dat)[[1]],
    nroutes=length(unique(all.dat$route)),
    nstrata=length(unique(all.dat$unit_code)),
    nstates=length(unique(all.dat$state_num)),
    nbcr=length(unique(all.dat$bcr_num2)), 
    nlcc= length(unique(all.dat$lcc_num2 )),
    ntime=length(unique(all.dat$surv_yr)) 
  )
}
names(dwrs) <- keeper.spp
# save(dwrs, all.dat, surv.dat, routes, uniq.time, file="data\\data.rdata")

###########
# CBC manip
###########
csp <- read_xlsx(".\\data\\BR-CBC_Circle_Species_Report_SQL_updated.xlsx")
ce1 <- read_xlsx(".\\data\\BR-CBC_Effort_Report_SQL_updated-1.xlsx")
ce2 <- read_xlsx(".\\data\\BR-CBC_Effort_Report_SQL_updated-2.xlsx")
cw <- read_xlsx(".\\data\\BR-CBC_Weather_Report_SQL.xlsx")
cdsys <- CRS(" +proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96
             +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs
             +ellps=GRS80 +towgs84=0,0,0")

csp$survID <- with(csp, paste(Abbrev, Count_yr))
# summarize CBC species data
csptab<- csp %>% pivot_wider(c(Abbrev, Count_yr), names_from=COM_NAME, values_from=how_many, 
                             values_fn=sum, values_fill=0 )
uniq.csp <- csp[!duplicated(csp$Abbrev), c(1:4,7)]
csptab <- merge(csptab, uniq.csp[, -5], by="Abbrev", no.dups=T)
# conversions
ce2$km <- ifelse(ce2$Description...7=="Miles", ce2$Distance*1.609, ce2$Distance)
ce2 <- ce2[ce2$Description...5=="car",]
# merge all the survey data
csurv1 <- merge(ce2, ce1, by=c("Abbrev", "Count_yr"))
csurv1 <- csurv1[!duplicated(paste(csurv1$Abbrev, csurv1$Count_yr)),] # duplicates are in Southern Cali, delete
csptab <- merge(csptab, csurv1, by=c("Abbrev", "Count_yr"))
# subset to surveys>= 1966 consistent with Sauer pubs
csptab <- csptab[ csptab$Count_yr>=105 & csptab$Count_yr<=120,]
names(csptab)[3:13] <- c("RTHA", "COHA", "AMKE", "WTKI", 
                         "NOHA", "PRFA", "RSHA", "GOEA", 
                         "RLHA", "BAEA", "FEHA")
csptab <- csptab[,-c(15,16,17,18,19,20,21,24,25,27,28,29,30,31,32)]
# Fill in nonsurveyed as NAs
# tabulate routes surveyed per strat and year
survs.cbc <- tapply(csptab$Abbrev, list(csptab$Abbrev, csptab$Count_yr) , length, default=NA)
# add NA surveys so all are included
mcbc <- melt(survs.cbc)
colnames(mcbc) <- c("Abbrev", "Count_yr", "count")
mcbc <- mcbc[is.na(mcbc$count),]
df.filler <- data.frame(Abbrev=mcbc$Abbrev, Count_yr=mcbc$Count_yr, 
                        RTHA=NA, COHA=NA, AMKE=NA, WTKI=NA,
                        NOHA=NA, PRFA=NA, RSHA=NA, GOEA=NA, 
                        RLHA=NA, BAEA=NA, FEHA=NA,
                        Name=NA,
                        Hours= mean(csptab$Hours, na.rm=T),
                        km= mean(csptab$km, na.rm=T),
                        Field_counters= mean(csptab$Field_counters, na.rm=T)
                        )
csptab <- rbind(csptab, df.filler)
# SPATIAL JOINS
circs <- cbind(circle_2013_map@data, circle_2013_map@coords)
csptab <- merge(csptab, circle_2013_map, by.x="Abbrev", by.y="unit_name")
# convert data frame to spatial 
csp.sp <- SpatialPointsDataFrame(coords = data.frame(csptab$coords.x1, csptab$coords.x2), 
                                 data = csptab,
                                 proj4string = crs(stratmod)) 
# merge with strata data
cbcst <- over(csp.sp, stratmod,  returnList = F)
csp.sp <- cbind(csp.sp, cbcst)
csp.sp<- csp.sp[!is.na(csp.sp$unit_name),]
# remove zeroes and NAs
csp.sp[is.na(csp.sp$km) , ] <- mean(csp.sp$km, na.rm=T)
csp.sp[csp.sp$km==0, ] <- mean(csp.sp$km, na.rm=T)
# subset to strata in wrs dat
csp.sp <- csp.sp[ csp.sp$unit_code %in% unique(all.dat$unit_code), ]
strat_sub<- stratmod[stratmod$unit_code %in% unique(all.dat$unit_code), ]

writeOGR(obj=csp.sp, dsn="C://Users//rolek.brian//Documents//Projects//northwest-wrs//data//shapefiles", 
         layer="cbc counts", driver="ESRI Shapefile")
writeOGR(obj=strat_sub, dsn="C://Users//rolek.brian//Documents//Projects//northwest-wrs//data//shapefiles", 
         layer="stratum_modified_clipped", driver="ESRI Shapefile")
###############
# prepare data for JAGS
################
# extract corresponding values from wrs data for indices
bcr_tab<- all.dat[!duplicated(all.dat$bcr_num2), c("bcr_num", "bcr_num2")]
lcc_tab<- all.dat[!duplicated(all.dat$lcc_num2), c("lcc_num", "lcc_num2")]
state_tab <- all.dat[!duplicated(all.dat$state_num), c("unit_name", "state_num")]
strata_tab <- all.dat[!duplicated(all.dat$strata_num), c("unit_code", "strata_num")]
# assign these values to cbcb data
change.values <- function(x,y){ # x is the ref table and y is the data values
                                z <- rep(NA, length(y))
                                for (i in 1:nrow(x)){
                                            ind <- which(y==x[i,1]) 
                                            z[ind] <- x[i,2] 
                                            } 
                                          return(z) } # i and function
csp.sp$bcr_num2 <- change.values(bcr_tab, csp.sp$bcr_num)
csp.sp$lcc_num2 <- change.values(lcc_tab, csp.sp$lcc_num)
csp.sp$state_num2 <- change.values(state_tab, csp.sp$unit_name)
csp.sp$strata_num2 <- change.values(strata_tab, csp.sp$unit_code )
# calculate how many detections occurred in each strata
keeper.cbc <- keeper.spp
keeper.cbc[4] <- "BAEA"
csp3 <- csp.sp[ !is.na(csp.sp$RTHA), ]
csp.sp <- csp.sp[csp.sp$Abbrev %in% csp3$Abbrev,]
csp.sp$route_num2 <- as.numeric(factor(csp.sp$Abbrev))
# calc first
ind3 <- table( csp3$Abbrev, csp3$Count_yr+1899) # from CBC metadata: count years start at 1 in the 1900-1901 count season.
which.fun <- function(x){  min(which(x>0)) }
first.surv.cbc <- apply(ind3, 1, which.fun)
# merge won't work so using a backasswards way of merging
csp.sp99 <- merge(csp.sp@data, first.surv.cbc, by.x="Abbrev", by.y=0, duplicateGeoms = TRUE)
csp.sp99$first.surv.yr <- csp.sp99$y+2003
csp.sp$first.surv.yr <- csp.sp99$first.surv.yr
frst.cbc <- csp.sp99[ (csp.sp99$Count_yr+1899)==csp.sp99$first.surv.yr , ]
csp.sp$time_num <- (csp.sp$Count_yr+1899)-2003 # 2004-2005 is yr 1
csp.sp$yr <- csp.sp$time_num+2003

ldat.cbc <- ldat.wrs <- dall <- dcbc <- list()

for (sp in 1:length(keeper.cbc)){
  # Restrict data to strata where species were detected at least once
  keeper.strats.cbc <- names(table(csp.sp$strata_num2[csp.sp@data[, keeper.cbc[sp]]>0]))
  keeper.strats.wrs <- names(table(all.dat$strata_num[all.dat[, keeper.spp[sp]]>0]))
  keeper.strats <- as.numeric(unique(c(keeper.strats.cbc, keeper.strats.wrs)))
  cbcsub <- csp.sp[csp.sp$strata_num2 %in% keeper.strats , ]
  frst.cbc.sub <- frst.cbc[frst.cbc$strata_num2 %in% keeper.strats, ]
  wrssub <- all.dat[all.dat$strata_num %in% keeper.strats , ]
  frst.wrs.sub <- frst.wrs[frst.wrs$route_num %in% unique(wrssub$route_num), ]
  #  order so nas are last then remove duplicates
  # trying to keep values if they exist then keep NAs if 
  # they are the only values
  wrssub<- wrssub[order(is.na(wrssub$total)),]
  wrssub <- wrssub[!duplicated(data.frame(wrssub$route_num , wrssub$time_num) ), ]
  cbcsub <- cbcsub[order(is.na(cbcsub@data[,keeper.cbc[sp]] )),]
  cbcsub <- cbcsub[!duplicated(data.frame(cbcsub@data$route_num2 , cbcsub@data$time_num )), ]
  # remove times before first
  cbcsub <- cbcsub[cbcsub$time_num >= (cbcsub$first.surv.yr-2003), ]
  wrssub <- wrssub[(wrssub$time_num+2003) >= wrssub$first.surv.yr, ]
  strats <- c(wrssub$strata_num, cbcsub$strata_num2)
  strats <- as.numeric(factor(strats))
  ldat.wrs[[sp]] <- wrssub
  ldat.cbc[[sp]] <- cbcsub
  
  # create indices for WRS summaries
  routes  <- wrssub[!duplicated(wrssub$route_num),]
  routes <- routes[order(routes$route_num), ]
  ind.bcr <- table(routes$bcr_num2)
  ind.lcc <- table(routes$lcc_num2)
  ind.strata <- table(routes$strata_num)
  ind.state <- table(routes$state_num)
  bcr <- array(NA, dim=c(length(ind.bcr), max(ind.bcr)) )
  lcc <- array(NA, dim=c(length(ind.lcc), max(ind.lcc)) )
  strata <- array(NA, dim=c(length(ind.strata), max(ind.strata)) )
  state <- array(NA, dim=c(length(ind.state), max(ind.state)) )
  
  for (z in 1:length(ind.bcr)){
    w <- which(routes$bcr_num2==names(ind.bcr)[z])
    bcr[z,1:length(w)] <- w
  }
    for (z in 1:length(ind.lcc)){
    w2 <- which(routes$lcc_num2==names(ind.lcc)[z])
    lcc[z,1:length(w2)] <- w2
  }
  for (zz in 1:length(ind.strata)){
    w <- which(routes$strata_num==names(ind.strata)[zz])
    strata[zz,1:length(w)] <- w
  }
  for (zzz in 1:length(ind.state)){
    w <- which(routes$state_num==names(ind.state)[zzz])
    state[zzz,1:length(w)] <- w
  }
  # create indices for CBC summaries
  routes2  <- cbcsub[!duplicated(cbcsub$route_num2),]
  routes2 <- routes2[order(routes2$route_num2), ]
  ind.bcr2 <- table(routes2$bcr_num2)
  ind.lcc2 <- table(routes2$lcc_num2)
  ind.strata2 <- table(routes2$strata_num2)
  ind.state2 <- table(routes2$state_num2)
  bcr2 <- array(NA, dim=c(length(ind.bcr2), max(ind.bcr2)) )
  lcc2 <- array(NA, dim=c(length(ind.lcc2), max(ind.lcc2)) )
  strata2 <- array(NA, dim=c(length(ind.strata2), max(ind.strata2)) )
  state2 <- array(NA, dim=c(length(ind.state2), max(ind.state2)) )
  
  for (z in 1:length(ind.bcr2)){
    w <- which(routes2$bcr_num2==names(ind.bcr2)[z])
    bcr2[z,1:length(w)] <- w
  }
  for (z in 1:length(ind.lcc2)){
    w2 <- which(routes2$lcc_num2==names(ind.lcc2)[z])
    lcc2[z,1:length(w2)] <- w2
  }
  for (zz in 1:length(ind.strata2)){
    w <- which(routes2$strata_num2==names(ind.strata2)[zz])
    strata2[zz,1:length(w)] <- w
  }
  for (zzz in 1:length(ind.state2)){
    w <- which(routes2$state_num2==names(ind.state2)[zzz])
    state2[zzz,1:length(w)] <- w
  } 
   
  dcbc[[sp]] <- list( 
    count2= cbcsub@data[ ,keeper.cbc[sp] ],
    dist2= cbcsub$km/100, # per 100 km
    time2= cbcsub$time_num,
    rt2= as.numeric(factor(cbcsub$route_num2)),
    frst2= frst.cbc.sub$y, # 1965-66 is yr 1
    mn1st2= frst.cbc.sub[, keeper.cbc[sp] ] , 
    str_ind2= frst.cbc.sub$strata_num2,
    nroutesInStrata2= table(frst.cbc.sub$strata_num2), 
    state_ind2=frst.cbc.sub$state_num2,
    nroutesInState2=table(frst.cbc.sub$state_num2),
    bcr_ind2=frst.cbc.sub$bcr_num2,
    nroutesInBCR2=table(frst.cbc.sub$bcr_num2),
    lcc_ind2=frst.cbc.sub$lcc_num2,
    nroutesInLCC2=table(frst.cbc.sub$lcc_num2), 
    #sp=all.dat$sp_std, 
    ncounts2=nrow(cbcsub),
    nroutes2=length(unique(cbcsub$route_num2)),
    nstrata2=length(unique(cbcsub$strata_num2)),
    nstates2=length(unique(cbcsub$state_num2)),
    nbcr2=length(unique(cbcsub$bcr_num2)), 
    nlcc2= length(unique(cbcsub$lcc_num2 )),
    ntime2=length(unique(cbcsub$time_num)) 
  )
  bothstrata <- c(frst.wrs.sub$strata_num, frst.cbc.sub$strata_num2)
  bothstrata2 <- as.numeric(factor(bothstrata))
  str <- data.frame(old=bothstrata[!duplicated(bothstrata)],
                    new = bothstrata2[!duplicated(bothstrata2)])
  strata_tab[,2+sp] <- str$new[match(strata_tab$strata_num, str$old)]
  colnames(strata_tab)[2+sp] <- keeper.spp[sp]
  wrs.strata <- strata_tab[,2+sp][match (frst.wrs.sub$strata_num, strata_tab$strata_num)]
  cbc.strata <- strata_tab[,2+sp][match (frst.cbc.sub$strata_num2, strata_tab$strata_num)]
  
  dall[[sp]] <- list(
    # WRS data
    count= wrssub[ ,keeper.spp[sp] ],
    dist= (wrssub$Miles*1.609)/100, # per 100 km
    time= wrssub$time_num,
    rt= as.numeric(factor(wrssub$route_num)),
    frst= frst.wrs.sub$first.surv.yr-2003 ,
    mn1st= frst.wrs.sub[ , keeper.spp[sp] ],
    strat= bothstrata2,
    frstall= c(frst.wrs.sub$first.surv.yr-2003, frst.cbc.sub$y),
    strat1= wrs.strata,
    str_ind= strata,
    nroutesInStrata= ind.strata, 
    state_ind= state,
    nroutesInState= ind.state,
    bcr_ind= bcr,
    nroutesInBCR= ind.bcr,
    lcc_ind= lcc,
    nroutesInLCC= ind.lcc, 
    sp= wrssub$sp_std, 
    ncounts= nrow(wrssub),
    nroutes= length(unique(wrssub$route)),
    nstrata= length(strata_tab[,2+sp][ !is.na(strata_tab[,2+sp]) ] ),
    nstates= length(unique(wrssub$state_num)),
    nbcr= length(unique(wrssub$bcr_num2)), 
    nlcc= length(unique(wrssub$lcc_num2 )),
    ntime= length(unique(wrssub$surv_yr)), 
    # CBC data
    count2= cbcsub@data[ ,keeper.cbc[sp] ],
    dist2= cbcsub$km/100, # per 100 km
    time2= cbcsub$time_num,
    rt2= as.numeric(factor(cbcsub$route_num2)),
    frst2= frst.cbc.sub$y,
    mn1st2= frst.cbc.sub[ , keeper.cbc[sp] ], # 2003 is yr 1
    strat2=cbc.strata, 
    str_ind2= strata2,
    nroutesInStrata2= ind.strata2, 
    state_ind2=state2,
    nroutesInState2=ind.state2,
    bcr_ind2=bcr2,
    nroutesInBCR2=ind.bcr2,
    lcc_ind2=lcc2,
    nroutesInLCC2=ind.lcc2,
    ncounts2=nrow(cbcsub),
    nroutes2=length(unique(cbcsub$route_num2)),
    nstrata2=length(unique(cbcsub$strata_num2)),
    nstates2=length(unique(cbcsub$state_num2)),
    nbcr2=length(unique(cbcsub$bcr_num2)), 
    nlcc2= length(unique(cbcsub$lcc_num2 )),
    ntime2=length(unique(cbcsub$time_num)) 
)
} # sp
names(dcbc) <- names(dall) <- keeper.cbc
save(dwrs=dwrs, lwrs=all.dat, surv.wrs=surv.dat, route.wrs=routes, uniq.time, 
     dcbc=dcbc, lcbc=csp.sp, surv.cbc=frst.cbc, 
     dall=dall, strata_tab,
     file="data\\data.rdata")
                   


#################
## PLOTS 1
################
# plot survey locations
plot(stratmod[stratmod@data$unit_code %in% 
                unique(all$unit_code),], 
     axes=T, lwd=2, col="gray90",  
     xlim=c(-2250000,-1000000), ylim=c(200000,1300000))
points(spdf, pch=3, col="blue", lwd=3)
sub <- (cbc$unit_code %in% unique(all$unit_code)) &
  (substr(cbc$unit_name,1,2) %in% c("ID","WA","OR", "CA"))
points(circle_2013_map[sub,], pch=1, col="purple", lwd=3)

sub2 <- (substr(cbc$unit_name,1,2)=="CA") & cbc$lcc_num == 5
sub2[is.na(sub2)] <- FALSE
points(circle_2013_map[sub2,], pch=1, col="purple", lwd=3) 

centroids <- centroid(stratmod)
cents <- cbind(centroids, stratmod@data)
cents <- cents[cents$unit_code %in% unique(all$unit_code),]
text(cents[,c(1,2)], 
     labels=cents$unit_code,
     cex=1.2)
writeOGR(obj=spdf, dsn="C://Users//rolek.brian//Documents//Projects//northwest-wrs//data//shapefiles", 
         layer="wrs counts spdf", driver="ESRI Shapefile")
# text(cbc$coords.x1[sub2], cbc$coords.x2[sub2],
#      labels=as.character(cbc$unit_code[sub2]), 
#      col="purple")

# plot one weird point where the
# state asssigned using centroid 
# is not the same as state of survey 
# points(-1721590  , 910451.4, col="red", cex=4)

# plot survey locations
plot(stratmod[stratmod@data$unit_code %in% 
                unique(all$unit_code),], 
     axes=T, lwd=2, col="gray90",  
     xlim=c(-2250000,-1000000), ylim=c(200000,1300000))
points(spdf, pch=3, col="blue", lwd=3)
sub <- (cbc$unit_code %in% unique(all$unit_code)) &
  (substr(cbc$unit_name,1,2) %in% c("ID","WA","OR", "CA"))
points(circle_2013_map[sub,], pch=1, col="purple", lwd=3)

centroids <- centroid(stratmod)
cents <- cbind(centroids, stratmod@data)
cents <- cents[cents$unit_code %in% unique(all$unit_code),]
text(cents[,c(1,2)], 
     labels=cents$unit_code,
     cex=1.2)

# text(cbc$coords.x1[sub2], cbc$coords.x2[sub2],
#      labels=as.character(cbc$unit_code[sub2]), 
#      col="purple")

# plot one weird point where the
# state asssigned using centroid 
# is not the same as state of survey 
# points(-1721590  , 910451.4, col="red", cex=4)


plot(strat_sub)
points(csp.sp, pch=3, cex=0.5, col="#7fc97f", lwd=0.5)
points(spdf, pch=4, cex=0.5, col="#beaed4", lwd=0.5 )

m <- as(csp.sp, "data.frame") 
ggplot(m, aes(Longitude, Latitude)) + 
  geom_point() + coord_equal() +
  geom_area(stratmod)

# make gif over time
size <- col <- list()
newcols <- data.frame(AMKEeff=rep(NA, nrow(csp.sp)), 
              RTHAeff=NA, RLHAeff=NA, BAEGeff=NA, NOHAeff=NA)
csp.sp <- cbind(csp.sp, newcols)

for (sp in 1:1){
# Calcs for color scales and size
csp.sp@data[, 37+sp ] <- csp.sp[,2+sp][[1]]/(csp.sp$km/1000)
counteff <- log(csp.sp[,37+sp][[1]]+0.001)
size[[sp]] <- (counteff-min(counteff))/ 
              (max(counteff)-min(counteff))
col[[sp]] <- cut(size[[sp]], breaks=10)#breaks=c(-0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.5, 1))
png_path <- file.path("C:\\Users\\rolek.brian\\Documents\\Projects\\northwest-wrs\\docs\\figs\\frame%03d.png")
png(png_path, width = 6, height = 6, units = "in", res = 200)
par(ask = FALSE)
kp <- unique(csp.sp$unit_code)[!is.na(unique(csp.sp$unit_code))]
for (ts in sort(unique(csp.sp$yr)) ) {
  name <- paste(c( "AMKE", "RTHA", "RLHA", "BAEG", "NOHA")[sp], 
                " ", as.character(ts))
  par(bty="n")
  plot(stratmod[stratmod@data$unit_code %in% kp ,], 
       axes=T, lwd=1, col=NA, main=name, 
       xaxt="n", yaxt="n", family="sans", font=1,
       xlim=c(-2250000,-1000000), ylim=c(250000,1350000))
  s <- size[[sp]][csp.sp$yr==ts]
  bird <- csp.sp[csp.sp$yr==ts,]
  points(bird, pch=16, cex=s*2,
         col=magma(10, direction=1)[as.numeric(col[[sp]])] )
  legend()
} # ts
dev.off()

png_files <- sprintf(png_path, 1:length(unique(csp.sp$yr)))
gif_file <- tempfile(fileext = ".gif")
gifski(png_files, gif_file, delay = 0.2, progress = FALSE)
unlink(png_files)
utils::browseURL(gif_file)
} # sp
names(size) <- names(col) <- c( "AMKE", "RTHA", "RLHA", "BAEG", "NOHA")


#######
# Trying sf to plot in ggplot
########
csf <- st_as_sf(csp.sp)
st <- st_as_sf(stratmod[stratmod@data$unit_code %in% kp ,])
ggplot() + 
  geom_sf(data=st, aes(fill=as.factor(bcr_num))) +
  scale_fill_discrete() +
  theme_bw() +
  geom_sf(data=csf, aes(col=American.Kestrel)) +
  xlab("Longitude") + ylab("Latitude") +
  transition_states(Species,
                  transition_length = 2,
                  state_length = 1)



geom_point(data = sites, aes(x = longitude, y = latitude), size = 4, 
           shape = 23, fill = "darkred") +
  coord_sf(xlim = c(-88, -78), ylim = c(24.5, 33), expand = FALSE)

#################
## PLOTS 2
################
m <- as(csp.sp, "data.frame") 
ggplot(m, aes(Longitude, Latitude)) + 
  geom_point() + coord_equal() +
  geom_area(stratmod)

# make gif over time


csp.sp$sizeRT <- csp.sp@data$`Red-tailed Hawk`/max(csp.sp@data$`Red-tailed Hawk`)*4

csp.sp$colRT <- cut(csp.sp$sizeRT, 10)
png_path <- file.path("C:\\Users\\rolek.brian\\Documents\\Projects\\northwest-wrs\\docs\\figs\\frame%03d.png")

png(png_path, width = 6, height = 6, units = "in", res = 100)
par(ask = FALSE)

for (ts in sort(unique(csp.sp$yr)) ) {
  name <- as.character(ts)
  plot(stratmod[stratmod@data$unit_code %in% 
                  unique(cbcst$unit_code),], 
       axes=T, lwd=1, col=NA, main=ts,
       xlim=c(-2250000,-1000000), ylim=c(200000,1300000))
  s <- csp.sp[csp.sp$yr==ts,]
  points(s, pch=16, cex=s$sizeRT,
         col=magma(10, direction=1)[as.numeric(s$colRT)] )
}
dev.off()

png_files <- sprintf(png_path, 1:length(unique(csp.sp$yr)))
gif_file <- tempfile(fileext = ".gif")
gifski(png_files, gif_file, delay = 0.2, progress = FALSE)
unlink(png_files)
utils::browseURL(gif_file)


points(spdf, pch=3, col="blue", lwd=3)
sub <- (cbc$unit_code %in% unique(all$unit_code)) &
  (substr(cbc$unit_name,1,2) %in% c("ID","WA","OR"))
points(circle_2013_map[sub,], pch=1, col="purple", lwd=3)

sub2 <- (substr(cbc$unit_name,1,2)=="CA") & cbc$lcc_num == 5
sub2[is.na(sub2)] <- FALSE
points(circle_2013_map[sub2,], pch=1, col="purple", lwd=3)