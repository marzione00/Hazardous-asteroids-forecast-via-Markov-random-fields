library(readxl)
library(PerformanceAnalytics)
library(corrplot)
library(ggcorrplot)
library(mgm)
library(gRim)
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
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(ggpubr)

#Data Loading

Asteroids2 <- read_excel("Dataset/Asteroids2.xlsx")
Asteroids2<-na.omit(Asteroids2)





#Asteroids2_h<-subset(Asteroids2,Hazardous==TRUE)
#Asteroids2_n<-subset(Asteroids2,Hazardous==FALSE)

#Asteroids_NH<- sample(3878,2800)

#Asteroids2_n<-Asteroids2[Asteroids_NH,]

#Asteroids_NH_frame<-Asteroids2[Asteroids_NH,]

#Asteroids_FINAL<-rbind(Asteroids2_h,Asteroids2_n)

#save(Asteroids_FINAL,file="Asteroids_FINAL.rda")

load("Asteroids_FINAL.rda")

Asteroids_FINAL[,2:21]<-scale(Asteroids_FINAL[,2:21])

Asteroids_FINAL_double<-as.data.frame(lapply(Asteroids_FINAL[,2:22], as.double))
Asteroids_FINAL_double["Hazardous"] <- as.factor(Asteroids_FINAL_double[,21])



#preliminary analysis


SS <- CGstats(Asteroids_FINAL_double)
SS
cc<-apply(SS$center,1,sd) /apply(SS$center,1,mean)
cc

M_cov <- cov(Asteroids_FINAL_double[,-21])
M_cor <- cor(Asteroids_FINAL_double[,-21])

chart.Correlation(Asteroids_FINAL_double , histogram=TRUE, pch=19)
testRes = cor.mtest(M_cor, conf.level = 0.95)
corrplot(M_cor , p.mat = testRes$p, method = 'color', diag = FALSE, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9, 
         insig = 'label_sig', pch.col = 'grey20', order = 'AOE')



res.famd <- FAMD(Asteroids_FINAL_double)
fviz_famd_var(res.famd, "quanti.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)
fviz_famd_var(res.famd, "quali.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)
fviz_famd_ind(res.famd, col.ind = "cos2", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)

heatmap(x=M_cor)


Hazardous_subset<-subset(Asteroids_FINAL,Hazardous==TRUE)
Not_Hazardous_subset<-subset(Asteroids_FINAL,Hazardous==FALSE)

ggdensity(Asteroids_FINAL,x="Eccentricity",rug = TRUE, color = "Hazardous",fill = "Hazardous",size=2)+theme_bw()
ggdensity(Asteroids_FINAL,x="Absolute_Magnitude",rug = TRUE, color = "Hazardous",fill = "Hazardous",size=2)+theme_bw()
ggdensity(Asteroids_FINAL,x="Min_Orbit_Intersection",rug = TRUE, color = "Hazardous",fill = "Hazardous",size=2)+theme_bw()


#Continuous variables analysis


msx<-cmod(~.^1.,data=Asteroids_FINAL_double[,-21])
msx2<-stepwise(msx,direction="forward",k=log(nrow(Asteroids_FINAL_double[,-21])),details=1)
plot.igraph(as(msx2,'igraph'),layout=layout.fruchterman.reingold(as(msx2,'igraph'), niter=300000), vertex.color="red")

#Mixed interaction analysis


fit_mgm <- mgm(data = Asteroids_FINAL[,-1],type = c(rep("g",20),"c"),levels = c(rep(1,20),2),k=2,lambdaSel = "CV",lambdaFolds= 10,ruleReg="AND",overparameterize = T)


qgraph::qgraph(fit_mgm$pairwise$wadj,
               layout = "spring", repulsion = 1.3,
               edge.color = fit_mgm$pairwise$edgecolor,
               nodeNames = colnames(Asteroids_FINAL[,-1] ),
               color = c(rep("lightblue",20),"purple"),
               legend.mode="style2", legend.cex=.4,
               vsize = 3.5, esize = 15)

#Forecast assessment

pred_obj <- predict(fit_mgm, Asteroids_FINAL[,-1])

Predicted_vs_real<-as.data.frame(pred_obj[["predicted"]][,21])

Predicted_vs_real<-as.data.frame(as.factor(Predicted_vs_real$`pred_obj[["predicted"]][, 21]`))

colnames(Predicted_vs_real)<-c("PRED")

levels(Predicted_vs_real$PRED) <- c('FALSE', 'TRUE')

conf<-data.frame(lapply(Asteroids_FINAL[,22], as.factor),Predicted_vs_real)





Confusion_matrix_ASTEROIDS<-confusionMatrix(conf$Hazardous,conf$PRED)

fourfoldplot(Confusion_matrix_ASTEROIDS$table,color = c("red","darkgreen"), main = "Mixed Interaction")

Confusion_matrix_ASTEROIDS

pred_numeric<-prediction(as.numeric(conf$PRED),as.numeric(conf$Hazardous))

roc_svm.perf <- performance(pred_numeric, measure = "tpr", x.measure = "fpr")

phi_Asteroids<-performance(pred_numeric, "phi")

phi_Asteroids@y.values[[1]][2]

plot(roc_svm.perf,cex.lab=1.5,yaxis.cex.axis=1.5,xaxis.cex.axis=1.5)
abline(a=0, b= 1)






