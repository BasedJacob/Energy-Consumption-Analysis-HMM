library(plyr)
library(dplyr)
library(lubridate)
library(stats)
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)

#Read in the large data set and scale it
mydata <- read.table("C:/Users/cps/Documents/R/CMPT318/Project/TermProjectData.txt",header=TRUE,sep=",")
numericData <- mydata[3:9]
scaledData <- scale(numericData, center = TRUE, scale = apply(numericData, 2, sd, na.rm = TRUE))

scaledDf <- data.frame(scaledData)
scaledDf$Date <- mydata$Date
scaledDf$Time <- mydata$Time
scaledDf$Date_time <- as.POSIXct(paste(scaledDf$Date, scaledDf$Time), format = "%d/%m/%Y %H:%M:%S")

#Selection of day and time window
days <- subset(scaledDf, wday(scaledDf$Date_time) == 1)
hours <- subset(days, hour(days$Date_time) >= 8 & hour(days$Date_time) < 11)

#Removing NA values
hours <- hours[!is.na(hours$Global_active_power),]
hours <- hours[!is.na(hours$Global_reactive_power),]
hours <- hours[!is.na(hours$Global_intensity),]
hours <- hours[!is.na(hours$Voltage),]
hours <- hours[!is.na(hours$Sub_metering_1),]
hours <- hours[!is.na(hours$Sub_metering_2),]
hours <- hours[!is.na(hours$Sub_metering_3),]

apply(hours, 2, function(x) any(is.na(x) | is.infinite(x)))

#PCA and 2 plots to illustrate it
pca <- prcomp(hours[1:7], scale=TRUE )
print(pca)
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
ggbiplot(pca, alpha = 0)
barplot(pca.var.per, main="Bar Plot", xlab="Principal Component", ylab="Percent Variation")

#Count number of columns present in each unique day
ndays <- aggregate(Global_intensity ~ Date, data = hours, FUN = length)

#Train and test split
set.seed(1)
l <- c(ndays$Global_intensity)
l1 <- l[1:108]
l2 <- l[109:154]
dt = sort(sample(nrow(hours), 19439))
train <- hours[dt,]
test <- hours[-dt,]


#Train and fit models from 4~20 states 
library("depmixS4")

modmulti_ntimes4 <- depmix(list(Global_intensity~1,Global_active_power~1), data=train, nstates=4,
                          family=list(gaussian(),gaussian()), ntimes=l1, na.allow=FALSE)
fmmulti_ntimes4 <- fit(modmulti_ntimes4)


modmulti_ntimes8 <- depmix(list(Global_intensity~1,Global_active_power~1), data=train, nstates=8,
                           family=list(gaussian(),gaussian()), ntimes=l1, na.allow=FALSE)
fmmulti_ntimes8 <- fit(modmulti_ntimes8)


modmulti_ntimes12 <- depmix(list(Global_intensity~1,Global_active_power~1), data=train, nstates=12,
                           family=list(gaussian(),gaussian()), ntimes=l1, na.allow=FALSE)
fmmulti_ntimes12 <- fit(modmulti_ntimes12)


modmulti_ntimes16 <- depmix(list(Global_intensity~1,Global_active_power~1), data=train, nstates=16,
                          family=list(gaussian(),gaussian()), ntimes=l1, na.allow=FALSE)
fmmulti_ntimes16 <- fit(modmulti_ntimes16)


modmulti_ntimes20 <- depmix(list(Global_intensity~1,Global_active_power~1), data=train, nstates=20,
                           family=list(gaussian(),gaussian()), ntimes=l1, na.allow=FALSE)
fmmulti_ntimes20 <- fit(modmulti_ntimes20)

#Test model based on the selected number of states
modmulti_test <- depmix(list(Global_intensity~1,Global_active_power~1), data=test, nstates=20,
                        family=list(gaussian(),gaussian()), ntimes=l2, na.allow=FALSE)
pars <- getpars(fmmulti_ntimes20)
copy <- setpars(modmulti_test, pars)
fb <- forwardbackward(copy, return.all=TRUE)

#Loglik and BIC plots
x <- seq(4,20,4)
plot(x,c(BIC(fmmulti_ntimes4),BIC(fmmulti_ntimes8),BIC(fmmulti_ntimes12),BIC(fmmulti_ntimes16),BIC(fmmulti_ntimes20)),ty="b", xlab='nstates', ylab='Performance',col='blue')
plot(x,c(logLik(fmmulti_ntimes4),logLik(fmmulti_ntimes8),logLik(fmmulti_ntimes12),logLik(fmmulti_ntimes16),logLik(fmmulti_ntimes20)),ty="b", xlab='nstates', ylab='Performance',col='red')



#Processing of anomalous dataset1
mydata1 <- read.table("C:/Users/cps/Documents/R/CMPT318/Project/dataWithAnomalies1.txt",header=TRUE,sep=",")
numericdata1 <- mydata1[3:9]
scaleddata1 <- scale(numericdata1, center = TRUE, scale = apply(numericdata1, 2, sd, na.rm = TRUE))

scaledDf1 <- data.frame(scaleddata1)
scaledDf1$Date <- mydata1$Date
scaledDf1$Time <- mydata1$Time
scaledDf1$Date_time <- as.POSIXct(paste(scaledDf1$Date, scaledDf1$Time), format = "%d/%m/%Y %H:%M:%S")

days1 <- subset(scaledDf1, wday(scaledDf1$Date_time) == 1)
hours1 <- subset(days1, hour(days1$Date_time) >= 8 & hour(days1$Date_time) < 11)

hours1 <- hours1[!is.na(hours1$Global_active_power),]
hours1 <- hours1[!is.na(hours1$Global_reactive_power),]
hours1 <- hours1[!is.na(hours1$Global_intensity),]
hours1 <- hours1[!is.na(hours1$Voltage),]
hours1 <- hours1[!is.na(hours1$Sub_metering_1),]
hours1 <- hours1[!is.na(hours1$Sub_metering_2),]
hours1 <- hours1[!is.na(hours1$Sub_metering_3),]
l3 <- c(rep(as.integer(180),48))

modmulti_d1 <- depmix(list(Global_intensity~1,Global_active_power~1), data=hours1, nstates=20,
                        family=list(gaussian(),gaussian()), ntimes=l3, na.allow=FALSE)
copy1 <- setpars(modmulti_d1, pars)
fb_d1 <- forwardbackward(copy1, return.all=TRUE)



#Processing of anomalous dataset2
mydata2 <- read.table("C:/Users/cps/Documents/R/CMPT318/Project/dataWithAnomalies2.txt",header=TRUE,sep=",")
numericdata2 <- mydata2[3:9]
scaleddata2 <- scale(numericdata2, center = TRUE, scale = apply(numericdata2, 2, sd, na.rm = TRUE))

scaledDf2 <- data.frame(scaleddata2)
scaledDf2$Date <- mydata2$Date
scaledDf2$Time <- mydata2$Time
scaledDf2$Date_time <- as.POSIXct(paste(scaledDf2$Date, scaledDf2$Time), format = "%d/%m/%Y %H:%M:%S")

days2 <- subset(scaledDf2, wday(scaledDf2$Date_time) == 1)
hours2 <- subset(days2, hour(days2$Date_time) >= 8 & hour(days2$Date_time) < 11)

hours2 <- hours2[!is.na(hours2$Global_active_power),]
hours2 <- hours2[!is.na(hours2$Global_reactive_power),]
hours2 <- hours2[!is.na(hours2$Global_intensity),]
hours2 <- hours2[!is.na(hours2$Voltage),]
hours2 <- hours2[!is.na(hours2$Sub_metering_1),]
hours2 <- hours2[!is.na(hours2$Sub_metering_2),]
hours2 <- hours2[!is.na(hours2$Sub_metering_3),]

modmulti_d2 <- depmix(list(Global_intensity~1,Global_active_power~1), data=hours2, nstates=20,
                      family=list(gaussian(),gaussian()), ntimes=l3, na.allow=FALSE)
copy2 <- setpars(modmulti_d2, pars)
fb_d2 <- forwardbackward(copy2, return.all=TRUE)

#Processing of anomalous dataset3
mydata3 <- read.table("C:/Users/cps/Documents/R/CMPT318/Project/dataWithAnomalies3.txt",header=TRUE,sep=",")
numericdata3 <- mydata3[3:9]
scaleddata3 <- scale(numericdata3, center = TRUE, scale = apply(numericdata3, 2, sd, na.rm = TRUE))

scaledDf3 <- data.frame(scaleddata3)
scaledDf3$Date <- mydata3$Date
scaledDf3$Time <- mydata3$Time
scaledDf3$Date_time <- as.POSIXct(paste(scaledDf3$Date, scaledDf3$Time), format = "%d/%m/%Y %H:%M:%S")

days3 <- subset(scaledDf3, wday(scaledDf3$Date_time) == 1)
hours3 <- subset(days3, hour(days3$Date_time) >= 8 & hour(days3$Date_time) < 11)

hours3 <- hours3[!is.na(hours3$Global_active_power),]
hours3 <- hours3[!is.na(hours3$Global_reactive_power),]
hours3 <- hours3[!is.na(hours3$Global_intensity),]
hours3 <- hours3[!is.na(hours3$Voltage),]
hours3 <- hours3[!is.na(hours3$Sub_metering_1),]
hours3 <- hours3[!is.na(hours3$Sub_metering_2),]
hours3 <- hours3[!is.na(hours3$Sub_metering_3),]
modmulti_d3 <- depmix(list(Global_intensity~1,Global_active_power~1), data=hours3, nstates=20,
                      family=list(gaussian(),gaussian()), ntimes=l3, na.allow=FALSE)
copy3 <- setpars(modmulti_d3, pars)
fb_d3 <- forwardbackward(copy3, return.all=TRUE)

