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
library(gRapHD)
library(glasso)
library(SIN)
library(gRain)
library(randomForest)
library(randomForestExplainer)
library(e1071)
library(MASS)

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

S.asteroids <- cov.wt(Asteroids_FINAL_double[,-21], method="ML")$cov
K.asteroids <- solve(S.asteroids)
round(100*K.asteroids)

PC.asteroids <-cov2pcor(S.asteroids)
round(PC.asteroids*100)
round(100*PC.asteroids)
threeshold <-0.1
Z<-abs(PC.asteroids)
Z[Z<threeshold] <- 0
diag(Z) <- 0
Z[Z>0] <- 1

g.thresh <- as(Z,"graphNEL")
thresh.asteroids <- cmod(edgeList(g.thresh ), data=Asteroids_FINAL_double[,-21])
thresh.asteroids
plot(as(thresh.asteroids,"igraph"))

psin.asteroids<-sinUG(S.asteroids,n=nrow(Asteroids_FINAL_double[,-21]))
plotUGpvalues(psin.asteroids)

Zp<-psin.asteroids
Zp[]<-0
Zp[psin.asteroids<0.05]<- 1
g.p05<-as(Zp,"graphNEL")

pval105.asteroids <-cmod(edgeList(g.p05),data=Asteroids_FINAL_double[,-21])
plot(as(pval105.asteroids,"igraph"))
gsin.asteroids <-as(getgraph(psin.asteroids,0.1),"graphNEL")
plot(as(gsin.asteroids,"igraph"))


#S.body<- cov.wt(Asteroids_FINAL_double[,-21])$cov
C.asteroids  <- cov2cor(S.asteroids )
res.lasso <- glasso(C.asteroids , rho=0.3)
AM <- res.lasso$wi!=0
diag(AM)<-FALSE
g.lasso <- as(AM,"graphNEL")
nodes(g.lasso)<-names(Asteroids_FINAL_double[,-21])
glasso.asteroids <- cmod(edgeList(g.lasso),data=Asteroids_FINAL_double[,-21])
plot(as(glasso.asteroids ,"igraph"))





commonedges.asteroids <-as(msx2,"graphNEL")
othermodels<-list(gsin.asteroids ,thresh.asteroids ,msx2)


othermodels<-lapply(othermodels,as,"graphNEL")

for(i in 1:length(othermodels)) {
  commonedges.asteroids  <- graph::intersection(commonedges.asteroids ,othermodels[[i]])
}

plot(as(commonedges.asteroids ,"igraph"))



#minForest model

bf<-minForest(Asteroids_FINAL_double,homog=TRUE,forbEdges=NULL,stat="LR")
plot(bf)
mbG<-stepw(model=bf,data=Asteroids_FINAL_double,exact=TRUE)
plot(mbG,cex.vert.label=0.8,numIter=6000,col.labels=c("red"),vert.hl=c(21),col.hl=c("blue"))


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

round(pred_obj[["probabilties"]][[21]],2)

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






#Random forest



Asteroids_FINAL_double_train<- sample(3552,2344)
Asteroids_FINAL_double_test<-Asteroids_FINAL_double[-Asteroids_FINAL_double_train,]


RF_perf_out<-tuneRF(Asteroids_FINAL_double_train[,-21],Asteroids_FINAL_double_train[,21], ntree=5000)
RF_perf_out<-data.frame(RF_perf_out)
rfor.planet <-randomForest(Hazardous~.,data=Asteroids_FINAL_double,subset=Asteroids_FINAL_double_train,localImp = TRUE,importance=TRUE,proximity=TRUE,ntry=4)
rfor.predict<-data.frame(predict(rfor.planet, Asteroids_FINAL_double_test, type = "class"))

var_imp_rforest<-data.frame(varImp(rfor.planet))
colnames(var_imp_rforest)<-c("Variable","Overall")
var_imp_rforest[,1]<-rownames(var_imp_rforest)
rownames(var_imp_rforest)<-seq(1:20)

ggplot(var_imp_rforest, aes(y=reorder(Variable,Overall),x=Overall,color="red")) + 
  geom_point() +
  geom_segment(aes(x=0,xend=Overall,yend=Variable)) +
  scale_color_discrete(name="Variable Group") +
  xlab("Overall importance") +
  ylab("Variable Name") + guides(color = FALSE, size = FALSE) + theme_bw()


plot(rfor.planet)
tree_plot<-data.frame(rfor.planet[["err.rate"]])
tree_plot[4]<-seq(1:500)
colnames(tree_plot)<-c("OOB","Not_habitable","Habitable","Trees")




ggplot() + geom_line(data = tree_plot, aes(x = Trees, y = OOB,color = "OOB") ) + 
  geom_line(data = tree_plot, aes(x = Trees, y = Not_habitable,color = "Not H") ) +
  geom_line(data = tree_plot, aes(x = Trees, y = Habitable,color = "H") )+labs(color = "Legend")+theme() + xlab('Trees') + ylab('Error')+theme_bw()


plot(rfor.planet)
legend("top", colnames(rfor.planet$err.rate), fill=1:ncol(rfor.planet$err.rate))
varImpPlot(rfor.planet)
proximityPlot(rfor.planet)
print(rfor.planet)

rfor.predict["Test"]<-as.factor(Asteroids_FINAL_double_test[,21])

colnames(rfor.predict)<-c("Predict","Test")


fourfoldplot(table(rfor.predict), color = c("red","darkgreen"),main = "Random Forest")

pred_for<-prediction(as.numeric(rfor.predict$Predict),as.numeric(rfor.predict$Test))

roc_for.perf <- performance(pred_for, measure = "tpr", x.measure = "fpr")

autoplot(roc_for.perf)+theme_bw()


#SVM


tune_svm_full.out<-tune(svm ,Hazardous~.,data=Asteroids_FINAL_double[Asteroids_FINAL_double_train,], type = 'C-classification',kernel="polynomial",
                        ranges =list(cost=(1:10),degree=(1:5)))
print(tune_svm_full.out)
perf_svm<-data.frame(tune_svm_full.out[["performances"]])

ggplot(perf_svm,aes(x=cost,y=degree, z=error))+geom_line(color="red",linetype="dashed")+geom_point(color="red")+theme_bw()

#X11(width=60, height=60)
#plot_ly(perf_svm[,1:3],x = ~cost, y = ~degree, z = ~error, type="scatter3d", mode="markers") 


svm.full <- svm(Hazardous~., data=Asteroids_FINAL_double[Asteroids_FINAL_double_train,],type = 'C-classification', kernel="polynomial",cost=10,degree=3,)

svm.predict_full<-data.frame(predict(svm.full,Asteroids_FINAL_double[-Asteroids_FINAL_double_train,],type = "class"))

svm.predict_full["T"]<-as.factor(Asteroids_FINAL_double_test[,21])

svm_fin_full<-data.frame(svm.predict_full,stringsAsFactors = TRUE)

colnames(svm_fin_full)<-c("Predict","Test")

caret::confusionMatrix(table(svm_fin_full))

fourfoldplot(table(svm_fin_full), color = c("red","darkgreen"),conf.level = 0, margin = 1, main = "SVM_FULL")

pred_svm_full<-prediction(as.numeric(svm_fin_full$Predict),as.numeric(svm_fin_full$Test))

roc_svm_full.perf <- performance(pred_svm_full, measure = "tpr", x.measure = "fpr")

phi_svm_full<-performance(pred_svm_full, "mi")

phi_svm_full@y.values

autoplot(roc_svm_full.perf)+theme_bw()



#QDA


qda.planet<- qda(Hazardous~., data=Asteroids_FINAL_double, subset=Asteroids_FINAL_double_train)

qda.prob<-data.frame(predict(qda.planet,Asteroids_FINAL_double[-Asteroids_FINAL_double_train,],type = "response"))
qda.prob<-qda.prob["class"]

qda_fin<-data.frame(qda.prob,stringsAsFactors = TRUE)
qda_fin["Test"]<-as.factor(Asteroids_FINAL_double[-Asteroids_FINAL_double_train,21])

colnames(qda_fin)<-c("Predict","Test")

caret::confusionMatrix(table(qda_fin))

fourfoldplot(table(qda_fin), color = c("red","darkgreen"),conf.level = 0, margin = 1, main = "QDA")

pred_qda<-prediction(as.numeric(qda_fin$Predict),as.numeric(qda_fin$Test))

roc_qda.perf <- performance(pred_qda, measure = "tpr", x.measure = "fpr")

phi_qda<-performance(pred_qda, "phi")

plot(phi_qda)

autoplot(roc_qda.perf)+theme_bw()


#Logistic



model <- glm(Hazardous~.,family=binomial(link='logit'),data=Asteroids_FINAL_double)

summary(model)

logistic.prob<-data.frame(predict(model ,Asteroids_FINAL_double[-Asteroids_FINAL_double_train,],type = "response"))
colnames(logistic.prob)<-c("P")
logistic.prob<- data.frame(ifelse(logistic.prob > 0.5, "1", "0"))
logistic.prob["T"]<-as.factor(Asteroids_FINAL_double[-Asteroids_FINAL_double_train,21])

colnames(logistic.prob)<-c("P","T")


fourfoldplot(table(logistic.prob), color = c("red","darkgreen"), main = "Logistic")

caret::confusionMatrix(table(logistic.prob))
