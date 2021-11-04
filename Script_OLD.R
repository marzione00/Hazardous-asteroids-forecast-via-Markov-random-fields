library(readxl)
library(PerformanceAnalytics)
library(corrplot)
library(ggcorrplot)
library(mgm)
library(bnlearn)
library(gRapHD)
library(qgraph)
library(Rgraphviz)
library(igraph)
library(gRbase)
library(RBGL)
library(caret)
library(pROC)
library(ROCR)
library(ggplot2)

Asteroids <- read_excel("Dataset/Asteroids_REF.xlsx")
Asteroids2 <- read_excel("Dataset/Asteroids2.xlsx")
whitelist <- read_excel("Whitelist.xlsx")

Asteroids_double<-as.data.frame(lapply(Asteroids[,2:11], as.double))

Asteroids_double["Hazardous"] <- as.factor(Asteroids_double[,10])

Asteroids2_double<-as.data.frame(lapply(Asteroids2[,2:23], as.double))

Asteroids2_double["Hazardous"] <- as.factor(Asteroids2_double[,22])

SS <- CGstats(Asteroids_double)

M_cov <- cov(Asteroids_double[,2:9])
M_cor <- cor(Asteroids_double[,2:9])

chart.Correlation(Asteroids[,2:13] , histogram=TRUE, pch=19)
#corrplot(Asteroids[,2:13])
heatmap(x=M_cor)
ggcorrplot(M_cor)

cc<-apply(SS$center,1,sd) /apply(SS$center,1,mean)


msx<-cmod(~.^1.,data=Asteroids_double[,-10])
msx2<-stepwise(msx,direction="forward",k=log(nrow(Asteroids_double[,-10])),details=1)
plot(msx2,"neato")

msx<-cmod(~.^.,data=Asteroids_double[,-10])
msx2<-stepwise(msx,details=1,"test")
plot(msx2,"neato")
qgraph(msx2)


Asteroids2_double<-na.omit(Asteroids2_double)

msx<-cmod(~.^1.,data=Asteroids2_double[,-c(18,22)])
msx2<-stepwise(msx,direction="forward",k=log(nrow(Asteroids2_double[,1:15])),details=1)
plot.igraph(as(msx2,'igraph'),layout=layout.fruchterman.reingold(as(msx2,'igraph'), niter=300000), vertex.color="red")

igraph.from.graphNEL(msx2)

fit_mgm <- mgm(data = Asteroids[1:4000,2:11],type = c(rep("g",9),"c"),levels = c(rep(1,9),2),k=2,lambdaSel = "CV",lambdaFolds= 10,ruleReg="AND",overparameterize = T)

Asteroids2<-na.omit(Asteroids2)

fit_mgm <- mgm(data = Asteroids2[1:3500,2:23],type = c(rep("g",21),"c"),levels = c(rep(1,21),2),k=2,lambdaSel = "CV",lambdaFolds= 10,ruleReg="AND",overparameterize = T)



qgraph(fit_mgm$pairwise$wadj,edge.color = fit_mgm $pairwise$edgecolor,layout = "spring",labels =  Asteroids$colnames)
pred_obj <- predict(fit_mgm, Asteroids[1:3500,2:11])
pred_obj[["errors"]]

pred_obj <- predict(fit_mgm, Asteroids2[,2:23])
pred_obj[["errors"]]

ciccio<-as.data.frame(pred_obj[["predicted"]][,22])

ciccio<-as.data.frame(as.factor(ciccio$`pred_obj[["predicted"]][, 22]`))

colnames(ciccio)<-c("PRED")

levels(ciccio$PRED) <- c('FALSE', 'TRUE')

ciccio

conf<-data.frame(lapply(Asteroids2[,23], as.factor),ciccio)







peppo<-confusionMatrix(conf$Hazardous,conf$PRED)

fourfoldplot(peppo$table)

peppo

pred_svm<-prediction(as.numeric(conf$PRED),as.numeric(conf$Hazardous))

roc_svm.perf <- performance(pred_svm, measure = "tpr", x.measure = "fpr")

phi_svm<-performance(pred_svm, "mi")

phi_svm@y.values


plot(roc_svm.perf,cex.lab=1.5,yaxis.cex.axis=1.5,xaxis.cex.axis=1.5)
abline(a=0, b= 1)

ggplot:autoplot(roc_svm.perf)+theme_bw()



Asteroids[2:13] <- lapply(Asteroids[2:13], as.numeric)


qgraph::qgraph(fit_mgm$pairwise$wadj,
               layout = "spring", repulsion = 1.3,
               edge.color = fit_mgm$pairwise$edgecolor,
               nodeNames = colnames(Asteroids[1:2267,2:11] ),
               color = c(rep("lightblue",9),"purple"),
               legend.mode="style2", legend.cex=.4,
               vsize = 3.5, esize = 15)


qgraph::qgraph(fit_mgm$pairwise$wadj,
               layout = "spring", repulsion = 1.3,
               edge.color = fit_mgm$pairwise$edgecolor,
               nodeNames = colnames(Asteroids2[2:23] ),
               color = c(rep("lightblue",21),"purple"),
               legend.mode="style2", legend.cex=.4,
               vsize = 3.5, esize = 15)

plot.igraph(as(msx2,'igraph'),layout=layout.fruchterman.reingold(as(msx2,'igraph'), niter=3000000), vertex.color="green")

wl = matrix(c("Mean_Motion", "Hazardous"), ncol = 2, byrow = TRUE)

dag = mmhc(Asteroids[1:4000,2:11],whitelist = wl)
plot(as(amat(dag),"graphNEL"),)
graphviz.plot(dag, shape = "ellipse")