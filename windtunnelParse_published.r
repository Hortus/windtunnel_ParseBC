#Wind tunnel parsing script with file looping

#Start with finalanalysis directory containing reps B and C with complete pitot data
#set working directory
setwd("D:/PhD/Videos/Wind Tunnel/batch2/parseBC")

#get all bending files
bendfiles <- list.files(pattern = "^[bend]")
#get all coeff files
coefffiles <- list.files(pattern = "^[coeff]")
#get all pitot files in pitot subdirectory
pitotfiles <- list.files(pattern = "^[pitot]")

#average the bend file. Get average value for every 24 rows
#Create a function to get bin means
BinMean <- function (vec, every, na.rm = FALSE) {
  n <- length(vec)
  x <- .colMeans(vec, every, n %/% every, na.rm)
  r <- n %% every
  if (r) x <- c(x, mean.default(vec[(n - r + 1):n], na.rm = na.rm))
  x
}

#read in stalker force values (expressed in N)
stalker <- read.csv("stalker.csv",header=FALSE)
#read in wind tunnel data csv 
windtunnel <- read.csv("windtunnelMain.csv",header=TRUE)
#dataframe for for loop output
outdfBC <- data.frame()


#loop through the files (68).  Read bend files into the bend dataframe, coeff files into the coeff dataframe, pitot files into the pitot dataframe, and stalker file into a vorce variable. 
for (j in 1:length(bendfiles)){
  #create the bend, coeff dataframes
  bend <- read.csv(bendfiles[j],header=FALSE)
  coeff <- read.csv(coefffiles[j],header=FALSE)
  #create the pitot dataframe
  pitot <- read.table(pitotfiles[j],header=TRUE,sep="\t")
  #force variable
  force <- stalker$V2[j]
  #pot number
  potnum <- stalker$V1[j]
  
  #Names of each column in bend file 1)frame 2)angles 3)f(x) at x halfheight point 4)x displacement at halfheight 5)x      displacement at quarterheight 6)x value at halfheight point for a given angle 7)y value at halfheight point for a given   angle 8)length of stem calculating by integrating along entire length of fitted curve 9)x value of the halfheight        point calculated by integrating to half the length of the curve fit 10)y value of the halfheight point calculated by     integrating to half the length of the curve fit.  11)projected area of the plant given the angle
  
  #bin each value
  #average angles
  angles <- BinMean(bend$V2, every = 24)
  #average f(x) at x halfheight point
  fx <- BinMean(bend$V3, every = 24)
  #x displacement at halfheight
  xdisp0.5 <- BinMean(bend$V4, every = 24)
  #x displacement at quarterheight
  xdisp0.25 <- BinMean(bend$V5, every = 24)
  #x value at halfheight point for a given angle
  xangle <- BinMean(bend$V6, every = 24)
  #y value at halfheight point for a given angle
  yangle <- BinMean(bend$V7, every = 24)
  #length of stem calculated by integrating along entire length of fitted curve
  lencalc <- BinMean(bend$V8, every = 24)
  #x value of the halfheight point calculated by integrating to half the length of the curve fit
  xcalc <- BinMean(bend$V9, every = 24)
  #y value of the halfheight point calculated by integrating to half the length of the curve fit
  ycalc <- BinMean(bend$V10, every = 24)
  #average area projections. Convert to sq meters from sq cm
  areas <- BinMean(bend$V11, every = 24)/10000
  #initial area value
  Ainit <- areas[1]
  
  #coeff value averaging.  Fields are as follows in coeff file 1)a 2)b
  #a coefficient 
  a <- BinMean(coeff$V1, every = 24)
  #b coefficient 
  b <- BinMean(coeff$V2, every = 24)
  
  #figure out a way to match the rows in both the angles and pitot data that is based on where they change.  Pitot data represents 0.66 seconds of measurement, 0.33 seconds of writing
  
  #get the lagged diff in values between elements in the angle and pitot vectors.
  #pitot differences will be large, consistent, and positive at the initial ramp up of the motor. When the length of positive indicies before a negative index for the sign difference in a vector is at least five, get that first negative index of the pitot vector.  Then go to the angle vector and find the index where the absolute value of the 5 preceeding indices was at least 2.5.  Align the pitot index with this position of the angle index.  
  
  signpitot <- sign(diff(pitot$X4492.ai2.mean))
  absangle <- abs(diff(angles))
  
  #loop through signpitot vector, starting at the 6th index
  pitot140index <- rep(NA,45) #get the index where motor stops increasing at 140rpm.
  n <- 7 #set starting index
  while (n <= 45){
    if ((signpitot[n] == -1) & (sum(signpitot[(n-6):(n-1)]) ==6)){  #if sign of the nth signpitot index is negative and the preceeding 6 are all positive
      #print(n)
      pitot140index[n] <- n
    } else {
      #align pitot tube index at 7
      pitot140index[n] <- 7
    }
    n <- n + 1 #add to starting index
  }
  
  #loop through abs angle vector, starting at the 3rd index
  angle140index <- NA #get the index where motor stops increasing at 140rpm. Allow vector to grow until truth is met
  #n <- 3 #set starting index
  n <- 15 #set starting index. Happens after the fit function adjusts and leaves are blown down tunnel
  while (n <= 45){
    if ((sum(absangle[n:(n+3)])) < 0.4){  #if the sum of the nth to the nth+3 abs angle differences are less than 0.4 degrees
      #print(n)
      angle140index[n] <- n
    } else {
      #assume the angle index change is stabalized by index 21 (frame ~500)
      angle140index[n] <- 21
    }
    n <- n + 1 #add to starting index
  }
  
  #Omit the na values
  pitot140index <- as.vector(na.omit(pitot140index))
  angle140index <- as.vector(na.omit(angle140index))
  
  #align the pitot dataframe with the angle vector at the first element of the pitot140 and angle140 indices
  
  if (pitot140index[1] > angle140index[2]){ #pitot index 1 is greater than angle index 1
    pitot2 <- pitot[((pitot140index[1]-angle140index[1]):(((pitot140index[1]-angle140index[1])+length(angles))-1)),] #subtract angle index from pitot index
  }else{ #pitot index 1 is less than angle index 1
    pitot2 <- pitot[((angle140index[1]-pitot140index[1]):(((angle140index[1]-pitot140index[1])+length(angles))-1)),] #subtact pitot index from angle index
  }
  
  #velocity calculation
  #atmospheric pressure by replicate (mb)
  #rep B 101-114 on 3/21, repB 115-139 on 3/22
  AtmPres <- 992.2
  #dew point. Consistent 1 degree dew point in the wind tunnel during the week
  DewP <- 1
  #Vapor pressure calculation (in mb)
  VapPres <- (6.1078*(10^((7.5*DewP)/(DewP+237.3)))) 
  #Pressure of dry air (in mb)
  PresDry <- AtmPres - VapPres
  #convert vapor pressure, dry air pressure to Pa
  VapPresPa <- VapPres*100
  PresDryPa <- PresDry*100
  #specific gas constants for dry air and water vapor
  Rair <- 287.05 # J/kg*K
  Rv <- 461.495 # J/kg*K
  #Density of mercury
  DensHg <- 13579 # kg/m3
  #calculate the density of the air at each frame, use this value along with the dynamic pressure (X4492.ai2.mean) to calculate the velocity at each frame
  Velocity <- c()
  DensAir <- c()
  for (i in 1:length(pitot2$Vertical)){
    TempK <- pitot2$Thermocouples.ai1.mean[i] + 273.15 #wind tunnel temp in Kelvin
    DP <- pitot2$X4492.ai2.mean[i]/10 #Dynamic pressure in mm Hg
    DensAir[i] <- ((PresDryPa/(Rair*TempK))+(VapPresPa/(Rv*TempK))) #kg/m3
    Velocity[i] <- sqrt((19.6*(DP/1000))*(DensHg/DensAir[i]))
  }
  
  #calculate the drag coefficient.  Base the velocity, density of the air on the index where angle first reached 50 degrees. 
  #first index for air density, velocity where angle index first reaches 50. Check indices beyond the 20th index in the angles vector
  z <- 20 #set starting index
  while (z <= 175){
    if (min(angles) < 50) {
      angleindex1 <- which(angles[20:length(angles)] <= 51 & angles[20:length(angles)] >=45)[1] #for pots where bending angle is reached
      angleindex <- angleindex1 + 20 #add 15 since searched vector length is shorter
    }else{
      angleindex1 <- which.min(angles[20:length(angles)]) #for pots where bending angle is not reached
      angleindex <- angleindex1 + 20 #add 15 since searched vector length is shorter
    }
    #add to index
    z <- z+1
  }
  #Estimate area of leaves at angleindex, and substract this value from the projected area at angleindex to get frontal area 
  Aleaf <- ((Ainit*100) / ((areas[angleindex])*100)^2) / 100
  Afront <- areas[angleindex] - Aleaf
  
  #use the values at this angle index to calculate drag coefficient 
  #cd <- (2*force)/(DensAir[angleindex]*(Velocity[angleindex]^2)*areas[angleindex])
  cd <- (2*force)/(DensAir[angleindex]*(Velocity[angleindex]^2)*Afront)
  
  # #option: keep only cd values where 50 degrees was reached (exculding first 15 indices of angles)
  # if (min(angles[15:length(angles)])<50){
  #   cd <- cd
  # }else{
  #   cd <- NA
  # }
  
  #now combine other variables into a dataframe, including velocity, model coefficients, displacement variables. 
  benddf <- data.frame(Velocity,angles,fx,xdisp0.5,xdisp0.25,xangle,yangle,lencalc,xcalc,ycalc,areas,a,b)
  
  #get values of a,b,angle,xdisp0.5,xdisp0.25 at the max velocity
  maxVelInd <- which.max(Velocity)
  max_a <- a[maxVelInd]
  max_b <- b[maxVelInd]
  max_angle <- angles[maxVelInd]
  max_xdisp0.5 <- xdisp0.5[maxVelInd]
  max_xdisp0.25 <- xdisp0.25[maxVelInd]
  max_area <- areas[maxVelInd]
  #get area,exponential coeff, drag force, and velocity where at index of drag coefficient estimation
  cd_area <- areas[angleindex]
  cd_velocity <- Velocity[angleindex]
  cd_b <- b[angleindex]
  force_d <- force
  #get the coefficient of lodging resistance at the point of drag coefficient estimation (angleindex)
  #estimate a value first. b is the known half height
  #half height vector for reps B,C
  halfheightBC <- c(101,	83,	94,	93,	100,	97,	96,	75,	110,	77,	109,	67,	99,	91,	99,	89,	84,	113,	94,	93,	116,	107,	86,	103,	95,	106,	108,	114,	93,	73,	83,	86,	73,	114,	87,	65,	71,	87,	87,	86,	84,	82,	71,	81,	97,	95,	82,	76,	81,	99,	102,	83,	88,	106,	95,	105,	99,	90,	92,	83,	81,	69,	86,	76,	75,	92,	100,	89)/2
  a <- xcalc[angleindex]*sqrt(1/(1-(ycalc[angleindex]^2/halfheightBC[j]^2)))
  CLr <- (a*force)/(halfheightBC[j]*xcalc[angleindex])
  #calculate the amount of angle recovery from max angle to end of video 
  recovery <- (angles[185]-angles[165])
  
  #append to outdf
  outdfBC <- rbind(outdfBC,c(potnum,cd,maxVelInd,max_a,max_b,max_angle,max_xdisp0.5,max_xdisp0.25,max_area,cd_area,cd_velocity,cd_b,force_d,CLr,recovery))
  
}

#change outdf colnames
colnames(outdfBC) <- c("potnum","cd","maxVelInd","max_a","max_b","max_angle","max_xdisp0.5","max_xdisp0.25","max_area","cd_area","cd_velocity","cd_b","force_d","CLr","recovery")


#Script for analyzing data from rep a, which does not contain pitot data. 
#set working directory
setwd("D:/PhD/Videos/Wind Tunnel/batch2/parseA")

#get all bending files
bendfiles <- list.files(pattern = "^[bend]")
#get all coeff files
coefffiles <- list.files(pattern = "^[coeff]")

#read in stalker force values (expressed in N)
stalkerA <- read.csv("stalker.csv",header=FALSE)

#dataframe for for loop output
outdfA <- data.frame()


#loop through the files (68).  Read bend files into the bend dataframe, coeff files into the coeff dataframe, pitot files into the pitot dataframe, and stalker file into a vorce variable. 
for (j in 1:length(bendfiles)){
  #create the bend, coeff dataframes
  bend <- read.csv(bendfiles[j],header=FALSE)
  coeff <- read.csv(coefffiles[j],header=FALSE)
  #create the pitot dataframe
  pitot <- read.table("pot114pitot.tsv",header=TRUE,sep="\t") #read from a representative pitot file.
  #force variable
  force <- stalkerA$V2[j]
  #pot number
  potnum <- stalkerA$V1[j]
  
  #Names of each column in bend file 1)frame 2)angles 3)f(x) at x halfheight point 4)x displacement at halfheight 5)x      displacement at quarterheight 6)x value at halfheight point for a given angle 7)y value at halfheight point for a given   angle 8)length of stem calculating by integrating along entire length of fitted curve 9)x value of the halfheight        point calculated by integrating to half the length of the curve fit 10)y value of the halfheight point calculated by     integrating to half the length of the curve fit.  11)projected area of the plant given the angle
  
  #bin each value
  #average angles
  angles <- BinMean(bend$V2, every = 24)
  #average f(x) at x halfheight point
  fx <- BinMean(bend$V3, every = 24)
  #x displacement at halfheight
  xdisp0.5 <- BinMean(bend$V4, every = 24)
  #x displacement at quarterheight
  xdisp0.25 <- BinMean(bend$V5, every = 24)
  #x value at halfheight point for a given angle
  xangle <- BinMean(bend$V6, every = 24)
  #y value at halfheight point for a given angle
  yangle <- BinMean(bend$V7, every = 24)
  #length of stem calculated by integrating along entire length of fitted curve
  lencalc <- BinMean(bend$V8, every = 24)
  #x value of the halfheight point calculated by integrating to half the length of the curve fit
  xcalc <- BinMean(bend$V9, every = 24)
  #y value of the halfheight point calculated by integrating to half the length of the curve fit
  ycalc <- BinMean(bend$V10, every = 24)
  #average area projections. Convert to sq meters from sq cm
  areas <- BinMean(bend$V11, every = 24)/10000
  #initial area value
  Ainit <- areas[1]
  
  #coeff value averaging.  Fields are as follows in coeff file 1)a 2)b
  #a coefficient 
  a <- BinMean(coeff$V1, every = 24)
  #b coefficient 
  b <- BinMean(coeff$V2, every = 24)
  
  #figure out a way to match the rows in both the angles and pitot data that is based on where they change.  Pitot data represents 0.66 seconds of measurement, 0.33 seconds of writing
  
  #get the lagged diff in values between elements in the angle and pitot vectors.
  #pitot differences will be large, consistent, and positive at the initial ramp up of the motor. When the length of positive indicies before a negative index for the sign difference in a vector is at least five, get that first negative index of the pitot vector.  Then go to the angle vector and find the index where the absolute value of the 5 preceeding indices was at least 2.5.  Align the pitot index with this position of the angle index.  
  
  signpitot <- sign(diff(pitot$X4492.ai2.mean))
  absangle <- abs(diff(angles))
  
  #loop through signpitot vector, starting at the 6th index
  pitot140index <- rep(NA,45) #get the index where motor stops increasing at 140rpm.
  n <- 7 #set starting index
  while (n <= 45){
    if ((signpitot[n] == -1) & (sum(signpitot[(n-6):(n-1)]) ==6)){  #if sign of the nth signpitot index is negative and the preceeding 6 are all positive
      #print(n)
      pitot140index[n] <- n
    } else {
      #align pitot tube index at 7
      pitot140index[n] <- 7
    }
    n <- n + 1 #add to starting index
  }
  
  #loop through abs angle vector, starting at the 3rd index
  angle140index <- NA #get the index where motor stops increasing at 140rpm. Allow vector to grow until truth is met
  #n <- 3 #set starting index
  n <- 15 #set starting index. Happens after the fit function adjusts and leaves are blown down tunnel
  while (n <= 45){
    if ((sum(absangle[n:(n+3)])) < 0.4){  #if the sum of the nth to the nth+3 abs angle differences are less than 0.4 degrees
      #print(n)
      angle140index[n] <- n
    } else {
      #assume the angle index change is stabalized by index 21 (frame ~500)
      angle140index[n] <- 21
    }
    n <- n + 1 #add to starting index
  }
  
  #Omit the na values
  pitot140index <- as.vector(na.omit(pitot140index))
  angle140index <- as.vector(na.omit(angle140index))
  
  #align the pitot dataframe with the angle vector at the first element of the pitot140 and angle140 indices
  
  if (pitot140index[1] > angle140index[2]){ #pitot index 1 is greater than angle index 1
    pitot2 <- pitot[((pitot140index[1]-angle140index[1]):(((pitot140index[1]-angle140index[1])+length(angles))-1)),] #subtract angle index from pitot index
  }else{ #pitot index 1 is less than angle index 1
    pitot2 <- pitot[((angle140index[1]-pitot140index[1]):(((angle140index[1]-pitot140index[1])+length(angles))-1)),] #subtact pitot index from angle index
  }
  
  #velocity calculation
  #atmospheric pressure by replicate (mb)
  #rep B 101-114 on 3/21, repB 115-139 on 3/22
  AtmPres <- 992.2
  #dew point. Consistent 1 degree dew point in the wind tunnel during the week
  DewP <- 1
  #Vapor pressure calculation (in mb)
  VapPres <- (6.1078*(10^((7.5*DewP)/(DewP+237.3)))) 
  #Pressure of dry air (in mb)
  PresDry <- AtmPres - VapPres
  #convert vapor pressure, dry air pressure to Pa
  VapPresPa <- VapPres*100
  PresDryPa <- PresDry*100
  #specific gas constants for dry air and water vapor
  Rair <- 287.05 # J/kg*K
  Rv <- 461.495 # J/kg*K
  #Density of mercury
  DensHg <- 13579 # kg/m3
  #calculate the density of the air at each frame, use this value along with the dynamic pressure (X4492.ai2.mean) to calculate the velocity at each frame
  Velocity <- c()
  DensAir <- c()
  for (i in 1:length(pitot2$Vertical)){
    TempK <- pitot2$Thermocouples.ai1.mean[i] + 273.15 #wind tunnel temp in Kelvin
    DP <- pitot2$X4492.ai2.mean[i]/10 #Dynamic pressure in mm Hg
    DensAir[i] <- ((PresDryPa/(Rair*TempK))+(VapPresPa/(Rv*TempK))) #kg/m3
    Velocity[i] <- sqrt((19.6*(DP/1000))*(DensHg/DensAir[i]))
  }
  
  #calculate the drag coefficient.  Base the velocity, density of the air on the index where angle first reached 50 degrees. 
  #first index for air density, velocity where angle index first reaches 50. Check indices beyond the 20th index in the angles vector
  z <- 20 #set starting index
  while (z <= 175){
    if (min(angles) < 50) {
      angleindex1 <- which(angles[20:length(angles)] <= 51 & angles[20:length(angles)] >=45)[1] #for pots where bending angle is reached
      angleindex <- angleindex1 + 20 #add 15 since searched vector length is shorter
    }else{
      angleindex1 <- which.min(angles[20:length(angles)]) #for pots where bending angle is not reached
      angleindex <- angleindex1 + 20 #add 15 since searched vector length is shorter
    }
    #add to index
    z <- z+1
  }
  
  #Estimate area of leaves at angleindex, and substract this value from the projected area at angleindex to get frontal area 
  Aleaf <- ((Ainit*100) / ((areas[angleindex])*100)^2) / 100
  Afront <- areas[angleindex] - Aleaf
  
  #use the values at this angle index to calculate drag coefficient 
  #cd <- (2*force)/(DensAir[angleindex]*(Velocity[angleindex]^2)*areas[angleindex])
  cd <- (2*force)/(DensAir[angleindex]*(Velocity[angleindex]^2)*Afront)
  
  # #optional write drag to file only if 50 degrees is reached. for first 15 indices of angles
  # if (min(angles[15:length(angles)])<50){
  #   cd <- cd
  # }else{
  #     cd <- NA
  # }

  #now combine other variables into a dataframe, including velocity, model coefficients, displacement variables. 
  benddf <- data.frame(Velocity,angles,fx,xdisp0.5,xdisp0.25,xangle,yangle,lencalc,xcalc,ycalc,areas,a,b)
  
  #get values of a,b,angle,xdisp0.5,xdisp0.25 at the max velocity
  maxVelInd <- which.max(Velocity)
  max_a <- a[maxVelInd]
  max_b <- b[maxVelInd]
  max_angle <- angles[maxVelInd]
  max_xdisp0.5 <- xdisp0.5[maxVelInd]
  max_xdisp0.25 <- xdisp0.25[maxVelInd]
  max_area <- areas[maxVelInd]
  #get area and velocity where at index of drag coefficient estimation
  cd_area <- areas[angleindex]
  cd_velocity <- Velocity[angleindex]
  cd_b <- b[angleindex]
  force_d <- force
  #get the coefficient of lodging resistance at the point of drag coefficient estimation (angleindex)
  #estimate a value first. b is the known half height
  #half height vector for rep A
  halfheightA <- c(86,	74,	84,	109,	95,	91,	93,	97,	99,	103,	85,	79,	83,	97,	72,	86,	114,	100,	84,	86,	66,	71,	72,	80,	91,	120,	93,	84)/2
  a <- xcalc[angleindex]*sqrt(1/(1-(ycalc[angleindex]^2/halfheightA[j]^2)))
  CLr <- (a*force)/(halfheightA[j]*xcalc[angleindex])
  #calculate the amount of angle recovery from angle bin 145 (12 m/s) to end of video 
  recovery <- (angles[185]-angles[165])
  
  
  #append to outdf
  outdfA <- rbind(outdfA,c(potnum,cd,maxVelInd,max_a,max_b,max_angle,max_xdisp0.5,max_xdisp0.25,max_area,cd_area,cd_velocity,cd_b,force_d,CLr,recovery))
  
}

colnames(outdfA) <- c("potnum","cd","maxVelInd","max_a","max_b","max_angle","max_xdisp0.5","max_xdisp0.25","max_area","cd_area","cd_velocity","cd_b","force_d","CLr","recovery")

#combine outdfA and outdfBC
outdf <- data.frame()
outdf <- rbind(outdf,outdfA,outdfBC)
#remove row 31 (pot 103) as the 2-row barley plant did not head and appeared diseased.  Drag coefficient was orders of magnitude larger than the others
outdf <- outdf[-31,]

#get merge biological and aerodynamic data
common <- intersect(windtunnel$pot,outdf$potnum)
finaldf <- data.frame()
for (i in 1:length(common)){
  commonTunnelind <- which(windtunnel$pot == common[i])
  commonOutdfInd <- which(outdf$potnum == common[i])
  finalrow <- cbind(windtunnel[commonTunnelind,],outdf[commonOutdfInd,])
  #append to final df
  finaldf <- rbind(finaldf,finalrow)
}



#write finaldf to file. 
write.csv(finaldf,file="Windtunnel_data.csv") 
  