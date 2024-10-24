## Winter 2022
## Project: Importation Risk
## Zahra Mohammadi

# Revised by Amy Hurford on October 18, 2024. General linear models were changed to
# generalized linear models (poisson). The original GLM code is commented out. The 
# mechanistic model is changed to be evaluated relative to a Poisson distribution.
# Rather than reporting R^2, the code was changed to report the likelihood test ratio
# statistic, since R^2 cannot be interpretted as the proportion of the variance
# explained for generalized linear models. A section was also added to consider travel
# volume with just the IATA data.


##---------- 1. Load libraries ------------
source("preliminaries.R")

##---------- 2. Read data ------------
# COVID cases in origins (CA excluded NL and US would consider as International)
source("getcases.pro.R") # output: Case, Cum and pop

#Seroprevalence Data and Correction factor calculated in 
source("getsero-CF.R") # output: CF

# daily travel volume calculated in Volume-focusOrigin.R
Volume <- read.csv("Data/TraVolData/TotalTraVolume.csv")
Volume$date <- as.Date(as.character(Volume$date), format = "%Y-%m-%d")
Volume <- Volume %>% rename(US = INT)

#break down daily volume to regular travel vs rotational workers (200 daily-rw)
rw_Volume <- Volume %>% mutate(AB = 114, ON= 30, BC= 8 , SK = 4, MB = 2, QC = 4, 
                               NB = 8, NS = 16, PEI = 4, TR = 10 , US = 0 )
#min(Volume$..) AB, Sk, PEI , NB and TR are less than this values needs to fix 
fixor <- c("AB","PEI", "SK", "TR", "NB")
for (i in 1:length(fixor)){
  pr <- fixor[i]
  for (j in 1:nrow(rw_Volume)){
    if (rw_Volume[j, pr] > Volume[j, pr] ){rw_Volume[j,pr] <- Volume[j, pr]}
  }
}

# Regular travelers
origins <- c("AB", "BC", "MB" ,"NB" ,"NS" , "ON","PEI", "QC","SK", "TR", "US"  )
rt_Volume <- Volume
for (i in 1:length(origins)){
  pr <- origins[i]
  rt_Volume[pr] <- Volume[pr] - rw_Volume[pr]
}

## Reported travel-related cases in NL
Travel_case <- read.csv("Data/Travel.Case.csv") 
Travel_case$date <- as.Date(Travel_case$date, format = "%Y-%m-%d") 

##-------------3. Estimate risk importation -----------------
# Note :## note a=0..13 need to read probsym and lambda =1:14 

# Effect correction factor(CF) from sero data on case data 
# to get the correct case data (for plotting)
CFCase <- function(df1, df2, origins){
  # inputs: df1 = correction , df2 = Case
  origins <- origins
  df <- df2
  for (i in 1:length(origins)){
    OR <- origins[i]
    for (j in 1:nrow(df)){
      df[j, OR] <- (df1[1, OR])*(df2[j, OR]) + df2[j, OR]
      }
    }
    return(as.data.frame(df))
}

origins <- c("AB", "BC","NB", "MB","NS","ON", "PEI","QC", "SK", "TR", "US")
Case2 <- CFCase(CF, Case, origins)

#Daily incidence proportion in departure (per date)
Dailyinc <- function(df1, pop, origins){
  # df1 = case2 , df2 = population
  df <- df1
  td <- 11       #delay between the exposed and report date
  origins <- origins
  for (i in 1:length(origins)){
    OR <- origins[i]
    pop2 <- pop$pop[pop$province == OR]
    for (j in 1:nrow(df)){
      td <- 11       #delay between the exposed and report date
      df[j, OR] <- df1[j+td, OR]/pop2}
  }
  return(df)
}
origins <- c("AB", "BC","NB", "MB","NS","ON", "PEI","QC", "SK", "TR", "US")
Dailyinfpre <- Dailyinc(Case2,pop, origins )
Dailyinfpre  <- Dailyinfpre  %>% filter(Dailyinfpre$date < as.Date("2021-06-01"))

## the rate of infected travelers including pre-test or not (changed value by 
## pre test date January 7, 2021 for international and May 15, 2021 for all)
T.INT <- function(t, a, df1, df2, pr){
  # df1 = Dailyinfpre , df2 = travel volume , pr = "US"
  rho <- 0.69 
  psi <- 0.75
  #PCR Test sensitivity - age infection 0-25
  tsens <- c(0, 0.05, 0.1, 0.55, 0.78, 0.77, 0.73, 0.68, 0.64, 0.59, 0.55, 
             0.49, 0.43,0.37, 0.31, 0.25, 0.22, 0.19, 0.16, 0.13, 0.1, 0.09,
             0.08, 0.07, 0.06, 0.05)
  #Probability develop symptom - age infection 0-16
  probsym = dgamma(seq(0,13), shape = 5.807, scale = 0.948) 
  lambda = cumsum(probsym)
  if (t < as.Date("2021-01-07")){
    x <- df2[df2$date == t, pr]*df1[df1$date == t, pr]*rho*psi*lambda[a]}
  else { 
    if (a > 2){x <- df2[df2$date == t, pr]*df1[df1$date == t, pr]*(1-tsens[a-2])}
    else {x <- 0}
      }
  return(x)
}


T.CA <- function(t, a, df1, df2, pr){
  # df1 = Dailyinfpre , df2 = travel volume (r or rw)
  rho <- 0.69 
  psi <- 0.75   
  #PCR Test sensitivity - age infection 0-25
  tsens <- c(0, 0.05, 0.1, 0.55, 0.78, 0.77, 0.73, 0.68, 0.64, 0.59, 0.55, 
             0.49, 0.43,0.37, 0.31, 0.25, 0.22, 0.19, 0.16, 0.13, 0.1, 0.09,
             0.08, 0.07, 0.06, 0.05)
  #Probability develop symptom - age infection 0-16
  probsym = dgamma(seq(0,13), shape = 5.807, scale = 0.948) 
  lambda = cumsum(probsym)
  ## note a=0..13 need to read probsym and lambda =1:14
  if (t < as.Date("2021-05-15")){
    x <- df2[df2$date == t, pr]*df1[df1$date == t, pr]*rho*psi*lambda[a]}
  else { 
    if (a > 2){x <- df2[df2$date == t, pr]*df1[df1$date == t, pr]*(1-tsens[a-2])}
    else {x <- 0}
  }
  return(x)
}

##################### International- regular travelers
R.INT.s = function(Dailyinfpre, rt_Volume, pr){
  df1 <- data.frame(date =seq(from=as.Date("2020-09-01"), to = as.Date("2021-05-31"), by = "days"))
  tr <- 3
  trep <- 1
  abar <- 1
  df1["report"] <- df1$date+ trep + abar + tr
  df1["R"] <- NA
  rho <- 0.69
  psi <- 0.75  
  PS <- 0.8
  #PCR Test sensitivity - age infection 0-25
  tsens <- c(0, 0.05, 0.1, 0.55, 0.78, 0.77, 0.73, 0.68, 0.64, 0.59, 0.55, 
            0.49, 0.43,0.37, 0.31, 0.25, 0.22, 0.19, 0.16, 0.13, 0.1, 0.09,
            0.08, 0.07, 0.06, 0.05)
  #Probability develop symptom - age infection 0-13
  probsym = dgamma(seq(0,13), shape = 5.807, scale = 0.948) 
  lambda = cumsum(probsym)
  for (i in 1:nrow(df1)){
    t <- df1$date[i]
    x <- 0
    for (a in 1:14){
      xa <- T.INT(t, a, Dailyinfpre, rt_Volume, pr)*tsens[a+abar+tr]
      x <- x + xa
      }
    df1$R[i] <- rho*PS*x
    }
  return(df1)
}

R.INT1 <- R.INT.s(Dailyinfpre, rt_Volume, "US") 


R.INT.e = function(rt_Volume, pr){
  df1 <- data.frame(date =seq(from=as.Date("2020-09-01"), to = as.Date("2021-05-31"), by = "days"))
  tr <- 3
  trep <- 1
  te <- 3
  df1["report"] <- df1$date + trep +  tr
  df1["E"] <- NA
  PA1 <- 0.01
  PA2 <- 0
  expo = c(as.Date("2020-09-27"), as.Date("2020-10-05"), as.Date("2020-10-23"),
           as.Date("2020-11-04"), as.Date("2020-11-23"), as.Date("2020-11-24"),
           as.Date("2020-12-04"), as.Date("2020-12-07"), as.Date("2020-12-11"),
           as.Date("2020-12-15"), as.Date("2020-12-17"), as.Date("2020-12-29"), 
           as.Date("2021-01-20"), as.Date("2021-02-15"), as.Date("2021-03-04"), 
           as.Date("2021-04-19"), as.Date("2021-05-02"), as.Date("2021-05-06"),
           as.Date("2021-05-11") )
  #PCR Test sensitivity - age infection 0-25
  tsens <- c(0, 0.05, 0.1, 0.55, 0.78, 0.77, 0.73, 0.68, 0.64, 0.59, 0.55, 
             0.49, 0.43,0.37, 0.31, 0.25, 0.22, 0.19, 0.16, 0.13, 0.1, 0.09,
             0.08, 0.07, 0.06, 0.05)
  #Probability develop symptom - age infection 0-13
  for (i in 1:nrow(df1)){
    t <- df1$date[i]
    if ( t  %in% expo) {
      df1$E[i] <- PA1*rt_Volume[rt_Volume$date == t, pr]*tsens[tr+te]
    }
    else {
      df1$E[i] <- PA2*rt_Volume[rt_Volume$date == t, pr]*tsens[tr+te]
    }
  }
  return(df1)
}

R.INT2 <- R.INT.e(rt_Volume, "US")
####################### Canadian Province -- regular travelers
R.CA.s = function(Dailyinfpre, rt_Volume, pr){
  df1 <- data.frame(date =seq(from=as.Date("2020-09-01"), to = as.Date("2021-05-31"), by = "days"))
  tr <- 3
  trep <- 1
  abar <- 1
  df1["report"] <- df1$date+ trep + abar + tr
  df1["R"] <- NA
  rho <- 0.69
  psi <- 0.75  
  PS <- 0.8
  #PCR Test sensitivity - age infection 0-25
  tsens <- c(0, 0.05, 0.1, 0.55, 0.78, 0.77, 0.73, 0.68, 0.64, 0.59, 0.55, 
             0.49, 0.43,0.37, 0.31, 0.25, 0.22, 0.19, 0.16, 0.13, 0.1, 0.09,
             0.08, 0.07, 0.06, 0.05)
  #Probability develop symptom - age infection 0-13
  probsym = dgamma(seq(0,13), shape = 5.807, scale = 0.948) 
  lambda = cumsum(probsym)
  for (i in 1:nrow(df1)){
    t <- df1$date[i]
    x <- 0
    for (a in 1:14){
      xa <- T.CA(t, a, Dailyinfpre, rt_Volume, pr)*tsens[a+abar+tr]
      x <- x + xa
    }
    df1$R[i] <-  rho*PS*x
  }
  return(df1)
}

# regular travelers is similar INT, so
Province <- function(Dailyinfpre, rt_Volume, origins){
  pro.df <- data.frame(date =seq(from=as.Date("2020-09-01"), to = as.Date("2021-05-31"), by = "days"))
  tr <- 3
  trep <- 1
  abar <- 1
  pro.df["report"] <- pro.df$date+ trep + abar + tr
  for (i in (1:length(origins))){
    pr <- origins[i]
    df <- R.CA.s(Dailyinfpre, rt_Volume, pr)  #dataframe output R.CA.s function
    pro.df[pr] <- df$R
  }
  return(pro.df)
}

origins <- c("ON", "QC","MB","NB","AB","NS","BC","SK","PEI","TR")
R.CA1 <- Province(Dailyinfpre, rt_Volume, origins )  

R.CA.e = function(rt_Volume, pr){
  df1 <- data.frame(date =seq(from=as.Date("2020-09-01"), to = as.Date("2021-05-31"), by = "days"))
  tr <- 3
  trep <- 1
  te <- 3
  df1["report"] <- df1$date + trep +  tr
  df1["E"] <- NA
  PA1 <- 0.01
  PA2 <- 0
  expo = c(as.Date("2020-09-27"), as.Date("2020-10-05"), as.Date("2020-10-23"),
           as.Date("2020-11-04"), as.Date("2020-11-23"), as.Date("2020-11-24"),
           as.Date("2020-12-04"), as.Date("2020-12-07"), as.Date("2020-12-11"),
           as.Date("2020-12-15"), as.Date("2020-12-17"), as.Date("2020-12-29"), 
           as.Date("2021-01-20"), as.Date("2021-02-15"), as.Date("2021-03-04"), 
           as.Date("2021-04-19"), as.Date("2021-05-02"), as.Date("2021-05-06"),
           as.Date("2021-05-11") )
  #PCR Test sensitivity - age infection 0-25
  tsens <- c(0, 0.05, 0.1, 0.55, 0.78, 0.77, 0.73, 0.68, 0.64, 0.59, 0.55, 
             0.49, 0.43,0.37, 0.31, 0.25, 0.22, 0.19, 0.16, 0.13, 0.1, 0.09,
             0.08, 0.07, 0.06, 0.05)
  #Probability develop symptom - age infection 0-13
  for (i in 1:nrow(df1)){
    t <- df1$date[i]
    if ( t  %in% expo) {
      df1$E[i] <- PA1*rt_Volume[rt_Volume$date == t, pr]*tsens[tr+te]
    }
    else {
      df1$E[i] <- PA2*rt_Volume[rt_Volume$date == t, pr]*tsens[tr+te]
    }
  }
  return(df1)
}

Province2 <- function(rt_Volume, origins){
  pro.df <- data.frame(date =seq(from=as.Date("2020-09-01"), to = as.Date("2021-05-31"), by = "days"))
  tr <- 3
  trep <- 1
  te <- 3
  pro.df["report"] <- pro.df$date + trep +  tr
  for (i in (1:length(origins))){
    pr <- origins[i]
    df <- R.CA.e(rt_Volume, pr)  #dataframe output R.CA.e function
    pro.df[pr] <- df$E
  }
  return(pro.df)
}

origins <- c("ON", "QC","MB","NB","AB","NS","BC","SK","PEI","TR")
R.CA2 <- Province2( rt_Volume, origins )  

######################## Canadian Province -- Rotational workers
RW.CA.s = function(j1, j2, j3, Dailyinfpre, rw_Volume, pr){
  df1 <- data.frame(date =seq(from=as.Date("2020-09-01"), to = as.Date("2021-05-31"), by = "days"))
  trep <- 1
  df1["report1"] <- df1$date+ trep + j1
  df1["S1"] <- NA
  df1["report2"] <- df1$date+ trep + j2
  df1["S2"] <- NA
  df1["report3"] <- df1$date+ trep + j3
  df1["S3"] <- NA
  #PCR Test sensitivity - age infection 0-25
  tsens <- c(0, 0.05, 0.1, 0.55, 0.78, 0.77, 0.73, 0.68, 0.64, 0.59, 0.55, 
             0.49, 0.43,0.37, 0.31, 0.25, 0.22, 0.19, 0.16, 0.13, 0.1, 0.09,
             0.08, 0.07, 0.06, 0.05)
  for (i in 1:nrow(df1)){
    t <- df1$date[i]
    x1 <- 0
    x2 <- 0
    x3 <- 0
    for (a in 1:14){
      xa1 <- T.CA(t, a, Dailyinfpre, rw_Volume, pr)*tsens[a+j1]
      x1 <- x1 + xa1
      xa2 <- T.CA(t, a, Dailyinfpre, rw_Volume, pr)*tsens[a+j2]*(1 - tsens[a+j1])
      x2 <- x2 + xa2
      xa3 <- T.CA(t, a, Dailyinfpre, rw_Volume, pr)*tsens[a+j3]*(1 - tsens[a+j1])*(1 - tsens[a+j2])
      x3 <- x3 + xa3
    }
    df1$S1[i] <-  x1
    df1$S2[i] <-  x2
    df1$S3[i] <-  x3
  }
  return(df1)
}

rw.AB <- RW.CA.s(3, 5, 7, Dailyinfpre, rw_Volume, "AB") # it is tested 
rw.ON <- RW.CA.s(3, 5, 7, Dailyinfpre, rw_Volume, "ON")
rw.PEI <- RW.CA.s(3, 5, 7, Dailyinfpre, rw_Volume, "PEI")
rw.QC <- RW.CA.s(3, 5, 7, Dailyinfpre, rw_Volume, "QC")
rw.SK <- RW.CA.s(3, 5, 7, Dailyinfpre, rw_Volume, "SK")
rw.BC <- RW.CA.s(3, 5, 7, Dailyinfpre, rw_Volume, "BC")
rw.MB <- RW.CA.s(3, 5, 7, Dailyinfpre, rw_Volume, "MB")
rw.NB <- RW.CA.s(3, 5, 7, Dailyinfpre, rw_Volume, "NB")
rw.NS <- RW.CA.s(3, 5, 7, Dailyinfpre, rw_Volume, "NS")
rw.TR <- RW.CA.s(3, 5, 7, Dailyinfpre, rw_Volume, "TR") # checked

RW.CA <- function(){
  df <- rw.AB 
  cols <- c("S1", "S2", "S3")
  df[, cols] <- rw.AB[, cols] + rw.ON[, cols] + rw.PEI[, cols] + rw.QC[, cols] +
                rw.SK[, cols] + rw.BC[, cols] + rw.MB[, cols] + rw.NB[, cols] + 
                rw.NS[, cols] + rw.TR[, cols]
  return(df)
}

RW.Pro <- RW.CA()
## merge three times in RW.Pro
RW.Pro2 <- function(){
  df1 <- select(RW.Pro, c(report1,S1)) %>% rename(report = report1)
  df2 <- select(RW.Pro, c(report2, S2)) %>% rename(report = report2)
  df3 <- select(RW.Pro, c(report3, S3)) %>% rename(report = report3)
  df1 <- df1 %>% left_join(df2, by = "report")
  df1 <- df1 %>% left_join(df3, by = "report")
  df1[is.na(df1)] <- 0
  df1["RW"] <- df1$S1 + df1$S2 + df1$S3
  df1 <- select(df1, c(report, RW))
  return(df1)
}

Imp.CA.rw <- RW.Pro2() 

#################aggregate
Imp.INT <- R.INT2 %>% left_join(R.INT1, by = "report")
Imp.INT <- select(Imp.INT, -c(date.x, date.y))
Imp.INT["RE.Int"] <- Imp.INT$E + Imp.INT$R
Imp.INT[is.na(Imp.INT)] <- 0

Imp.CA <- function(R.CA1, R.CA2){
  R.CA <- R.CA2 %>% left_join(R.CA1, by ="report")
  R.CA <- select(R.CA , -c(date.x, date.y))
  R.CA["ON"] <- R.CA$ON.x + R.CA$ON.y
  R.CA["QC"] <- R.CA$QC.x + R.CA$QC.y
  R.CA["MB"] <- R.CA$MB.x + R.CA$MB.y
  R.CA["NB"] <- R.CA$NB.x + R.CA$NB.y
  R.CA["AB"] <- R.CA$AB.x + R.CA$AB.y
  R.CA["NS"] <- R.CA$NS.x + R.CA$NS.y
  R.CA["BC"] <- R.CA$BC.x + R.CA$BC.y
  R.CA["SK"] <- R.CA$SK.x + R.CA$SK.y
  R.CA["PEI"] <- R.CA$PEI.x + R.CA$PEI.y
  R.CA["TR"] <- R.CA$TR.x + R.CA$TR.y
  R.CA <- select(R.CA, c("report", "ON", "QC","MB","NB","AB","NS","BC","SK","PEI","TR"))
  return(R.CA)
}

Imp.CA.r <- Imp.CA(R.CA1, R.CA2)  
Imp.CA.r[is.na(Imp.CA.r)] <- 0
Imp.CA.r["RE.Ca"] <- rowSums(Imp.CA.r[, c("ON", "PEI", "QC", "SK", "AB", "BC", "MB", "NB", "NS", "TR")])
# merger regular and rotational workers together and sum
Imp.CA <- Imp.CA.r %>% left_join(Imp.CA.rw, by = "report")
Imp.CA["REW.Ca"] <- rowSums(Imp.CA[, c("RE.Ca","RW")])

## one datafram for plotting
IMP <- Imp.CA %>% left_join(select(Imp.INT , c(report, RE.Int)) , by = "report")
IMP["Total"] <- rowSums(IMP[, c("REW.Ca","RE.Int")])
## Remove September5, since I just have rw not r 
IMP <- IMP %>% filter(IMP$report >= as.Date("2020-09-06") & IMP$report < as.Date("2021-06-01"))

##--------------- 5. competitive models (daily) -----------------------------------------
############# Models for international
df.model6 <- select( Travel_case, c(date, Tot.Int)) %>% left_join(Volume , by ="date")
df.model6 <- select(df.model6, c(date, Tot.Int, US)) %>% rename( vol = US, Rep.Int = Tot.Int )
df.model6 <- df.model6 %>% left_join(Dailyinfpre, by ="date")
df.model6 <- select(df.model6, c(date, Rep.Int, vol, US) ) %>% rename( infpre = US)
df.model6["volinf"] <- df.model6$vol*df.model6$infpre

#travel volume (as GLM)
#TVInt <- lm(Rep.Int ~ 0 + vol, data = df.model6)
#summary(TVInt)

TVInt <- glm(Rep.Int ~ 0 + vol, data = df.model6,family = poisson(link = "log"))
summary(TVInt)


#infection (as GLM is commented out)
#InfInt <- lm(Rep.Int ~ 0 + infpre, data = df.model6)
#summary(InfInt)

# As Poisson regression
InfInt <- glm(Rep.Int ~ 0 + infpre, data = df.model6, family = poisson(link = "log"))
summary(InfInt)


#both travel volume and infection (as GLM commented out)
#TVInfInt <- lm(Rep.Int ~ 0 + volinf , data = df.model6)
#summary(TVInfInt)

#both travel volume and infection
TVInfInt <- glm(Rep.Int ~ 0 + volinf , data = df.model6, family = poisson(link="log"))
summary(TVInfInt)


mod <- list(TVInt, InfInt, TVInfInt)
mod.name <-  c("TravelVolume-Int", "Infection.Prevalence-Int", "TravelVolume*Infection.Prevalence-Int")
AIC.INT = aictab(cand.set = mod, modnames = mod.name , c.hat = 1) 
LR.INT = -2*(-TVInt$null.deviance - c(logLik(TVInt),logLik(TVInfInt),logLik(InfInt)))
AIC.INT = data.frame(AIC.INT, LR = LR.INT)

######### Models for Domestic
df.model7 <- select( Travel_case, c(date, Tot.Dom)) %>% left_join(Volume , by ="date") 
df.model7 <- select(df.model7, -c(US)) %>% rename(Rep.Dom = Tot.Dom)
df.model7["vol.CA"] <- rowSums(df.model7[, c("ON", "PEI", "QC", "SK", "AB", "BC", "MB", "NB", "NS", "TR")])
#travel volume provinces (commented out as a GLM)
#TVPr <- lm(Rep.Dom ~ 0+ AB + BC + MB + NB + NS + ON + PEI + QC + SK + TR, data = df.model7 ) 
#summary(TVPr)
TVPr <- glm(Rep.Dom ~ 0+ AB + BC + MB + NB + NS + ON + PEI + QC + SK + TR, data = df.model7, family=poisson(link = "log"))
summary(TVPr)


#travel volume Canada (commented out as a GLM)
#TVCA <- lm( Rep.Dom ~ 0 + vol.CA, data= df.model7)
#summary(TVCA)
TVCA=glm( Rep.Dom ~ 0 + vol.CA, data= df.model7, family = poisson(link= "log"))
summary(TVCA)

df.model8 <- select( Travel_case, c(date, Tot.Dom)) %>% left_join( Dailyinfpre, by ="date") 
df.model8 <- select(df.model8, -c(US)) %>% rename(Rep.Dom = Tot.Dom)
df.model8["inf.CA"] <- rowSums(df.model8[, c("ON", "PEI", "QC", "SK", "AB", "BC", "MB", "NB", "NS", "TR")])
#infection provinces
#InfPr <- lm(Rep.Dom ~  0 + AB + BC + MB + NB + NS + ON + PEI + QC + SK + TR, data = df.model8)
#summary(InfPr)
InfPr <- glm(Rep.Dom ~  0 + AB + BC + MB + NB + NS + ON + PEI + QC + SK + TR, data = df.model8, family = poisson(link = "log"))
summary(InfPr)

#infection Canada
#InfCA <- lm(Rep.Dom  ~ 0+ inf.CA , data = df.model8 )
#summary(InfCA)
InfCA <- glm(Rep.Dom  ~ 0+ inf.CA , data = df.model8, family = poisson(link = "log"))
summary(InfCA)

##both travel volume*infection
multivolinf <- function(df1, df2, origins){
  # inputs: df1 = volume , df2 = infe
  origins <- origins
  df <- df1
  for (i in 1:length(origins)){
    OR <- origins[i]
    for (j in 1:nrow(df)){
      df[j, OR] <- df1[j, OR]*df2[j, OR]
    }
  }
  return(as.data.frame(df))
}
origins <- c("AB", "BC","NB", "MB","NS","ON", "PEI","QC", "SK", "TR")
df.model9 <- multivolinf(df.model7, df.model8, origins)
df.model9["volinfCA"] <- (df.model7$vol.CA)*(df.model8$inf.CA)
df.model9 <- select(df.model9, -c(vol.CA))
#province  
#TVInfPr <- lm(Rep.Dom  ~ 0 + AB + BC + MB + NB + NS + ON + PEI + QC + SK + TR ,  data = df.model9 )
#summary(TVInfPr)
TVInfPr=glm(Rep.Dom  ~ 0 + AB + BC + MB + NB + NS + ON + PEI + QC + SK + TR ,  data = df.model9 ,family = poisson(link = "log"))
summary(TVInfPr)


#Canada
#TVInfCA <- lm(Rep.Dom  ~ 0 + volinfCA , data = df.model9 )
#summary(TVInfCA)
TVInfCA <- glm(Rep.Dom  ~ 0 + volinfCA , data = df.model9 , family = poisson(link = "log"))
summary(TVInfCA)

models <- list(TVPr, TVCA, InfPr, InfCA, TVInfPr, TVInfCA )
model.names <-  c("Travel.Volume-Provinces", "Travel.Volume-CA", "Infection.Prevalence-Provinces", 
                  "Infection.Prevalence-CA", "TravelVolume*Infection.Prevalence-Provinces",
                  "TravelVolume*Infection.Prevalence-CA")
LR.CA = -2*(-TVPr$null.deviance - c(logLik(TVPr),logLik(TVCA),logLik(InfPr), logLik(InfCA), logLik(TVInfPr),logLik(InfCA)))
AIC.CA = aictab(cand.set = models, modnames = model.names, c.hat = 1)
AIC.CA = data.frame(AIC.CA,LR=LR.CA[c(5,3,1,6,2,4)])



################ model for epidemiology model output
df.modelfull <- select(IMP, c("report", "REW.Ca", "RE.Int")) %>% rename(date = report)
df.modelfull <- df.modelfull %>% left_join(select(Travel_case, c(date, Tot.Dom, Tot.Int)), by = "date") 
df.modelfull <- df.modelfull %>% rename(Est.Dom = REW.Ca, Est.Int = RE.Int ,
                                        Rep.Dom = Tot.Dom, Rep.Int = Tot.Int)
# commented out due to normal distribution
#set.seed(1066); 
#normalF <- function(sigma) {
#  x <- df.modelfull$Est.Dom  #Est.Int
#  y <- df.modelfull$Rep.Dom  #Rep.Int
#  sum (- 0.5*log(2*pi) - 0.5*log(sigma) - 0.5*(x - y)^2/sigma )
#}

#MLE = optim(c(0.1),               
#            fn = normalF,         # function to maximize
#            method = "L-BFGS-B",  # this method lets set lower bounds
#            lower = -Inf,         
#            control = list(fnscale = -1), # maximize
#            hessian = T                   # calculate Hessian matricce 
#)


##MLE
#MLE$value

##nll
#2*1 - 2*(MLE$value)

#R2
#R2(df.modelfull$Est.Int, df.modelfull$Rep.Int)

#RSQUARE = function(actual, preds){
#  sse = crossprod(actual - preds)
#  sst = crossprod(actual)
#  rsq <- 1-sse/sst
#  return(rsq)
#}
#RSQUARE(df.modelfull$Rep.Dom, df.modelfull$Est.Dom)

# No fitted parameters, calculated loglikelihood for mechanistic model
LL.CA=sum(dpois(df.modelfull$Rep.Dom, df.modelfull$Est.Dom, log = TRUE))
LL.INT=sum(dpois(df.modelfull$Rep.Int, df.modelfull$Est.Int, log = TRUE))


mech.AICc.CA = -2*LL.CA
mech.AICc.INT = -2*LL.INT

mech.mod.CA = c(Modnames = "mechanistic", K=0, AICc = mech.AICc.CA, Delta_AICc = mech.AICc.CA- as.numeric(AIC.CA$AICc[1]), ModelLik="", AICcWt = "", LL = LL.CA, Cum.Wt = "",LR="")
mech.mod.INT = c(Modnames = "mechanistic", K=0, AICc = mech.AICc.INT, Delta_AICc = mech.AICc.INT- as.numeric(AIC.INT$AICc[1]), ModelLik="", AICcWt = "", LL = LL.INT, Cum.Wt = "",LR="")

AIC.CA=rbind(AIC.CA, mech.mod.CA)
AIC.INT=rbind(AIC.INT, mech.mod.INT)

######### Save the best model output for both international and domestic
bestmodel.Int <- TVInt$fitted.values
bestmodel.CA <- TVInfPr$fitted.values
# make one dataframe including epi model, best model and reported for plotting 
# these are 273 values start from Sep 1 whereas predicted Imp start from Sep 6, 
# do adjust to make dataframe
dfplt <- df.modelfull
dfplt["bestInt"] <- bestmodel.Int[6:273]   
dfplt["bestDom"] <- bestmodel.CA[6:273]  
#there are couple negative in Dom
dfplt$bestDom[dfplt$bestDom < 0 ] <- 0

print(AIC.CA)
print(AIC.INT)
