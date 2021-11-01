library(readxl)
library(gRim)
library(PerformanceAnalytics)
library(corrplot)
library(ggcorrplot)

Asteroids <- read_excel("Dataset/Asteroids_REF.xlsx")
M_cov <- cov(Asteroids[,2:14])
M_cor <- cor(Asteroids[,2:14])

chart.Correlation(Asteroids[,2:13] , histogram=TRUE, pch=19)
corrplot(Asteroids[,2:13])
heatmap(x=M_cor)
ggcorrplot(M_cor)