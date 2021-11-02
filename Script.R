library(readxl)
library(gRim)
library(PerformanceAnalytics)
library(corrplot)
library(ggcorrplot)
library(mgm)
library(bnlearn)
library(gRapHD)
library(qgraph)

Asteroids <- read_excel("Dataset/Asteroids_REF.xlsx")
Asteroids[sapply(Asteroids, is.character)]  <- lapply(Asteroids[sapply(Asteroids, is.character)], as.factor)
Asteroids[2:8] <- lapply(Asteroids[2:8], as.numeric)

SS <- CGstats(Asteroids,varnames =c("Min_DIA","Max_DIA","Speed","MISS_DIST","Eccentricity","Semi_Major_Axis","Inclination","Orbital_Period","Perihelion_Distance","Aphelion_Dist","Perihelion_Time","Mean_motion","Hazardous"))


SS <- CGstats(Asteroids,varnames =c("Min_DIA","Speed","Eccentricity","Inclination","Orbital_Period","Hazardous"))

bF <- minForest(Asteroids)

M_cov <- cov(Asteroids[,2:14])
M_cor <- cor(Asteroids[,2:14])

chart.Correlation(Asteroids[,2:13] , histogram=TRUE, pch=19)
#corrplot(Asteroids[,2:13])
heatmap(x=M_cor)
ggcorrplot(M_cor)

cc<-apply(SS$center,1,sd) /apply(SS$center,1,mean)


msx<-dmod(~.^.,data=Asteroids[1:30,2:7])
msx2<-stepwise(msx,crit="test", alpha=0.15,details = 2)
plot(msx2)


fit_mgm <- mgm(data = Asteroids[1:4000,2:11],type = c(rep("g",9),"c"),levels = c(rep(1,9),2),k=2,lambdaSel = "CV",lambdaFolds= 10,ruleReg="AND",overparameterize = T)

dfnum = bnlearn.df2onehot(Asteroids[2:14])
plot(res)



qgraph(fit_mgm$pairwise$wadj,edge.color = fit_mgm $pairwise$edgecolor,layout = "spring",labels =  Asteroids$colnames)
pred_obj <- predict(fit_mgm, Asteroids[1:2267,2:11])
pred_obj[["errors"]]

Asteroids[2:13] <- lapply(Asteroids[2:13], as.numeric)


qgraph::qgraph(fit_mgm$pairwise$wadj,
               layout = "spring", repulsion = 1.3,
               edge.color = fit_mgm$pairwise$edgecolor,
               nodeNames = colnames(Asteroids[1:2267,2:11] ),
               color = c("purple",rep("lightblue",13)),
               legend.mode="style2", legend.cex=.4,
               vsize = 3.5, esize = 15)

dag = hc(Asteroids[2:14])
fit = bn.fit(dag, Asteroids[2:14])