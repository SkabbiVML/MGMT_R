Metsum_ROI_long <- read.csv("Data/Metsum_ROI_long.csv", header = T, row.names = 1)
colnames(Metsum_ROI_long) <- c(-7:-1, 1:98,c("+1","+2", "+3","+4", "+5", "+6", "+7", "+8", "+9", "+10"))
Metsum_ROI_long_ISLAND <- Metsum_ROI_long[c(2:nrow(Metsum_ROI_long)),c(8:105)]

load("Data/DotplotData.RData")

library(mixOmics)
library(RColorBrewer)
library(dplyr)

MGMT_RunSum <- as_tibble(MGMT_RunSum)

Training <- MGMT_RunSum %>% 
  filter(Series != "DenStem") %>% 
  filter(Pyro_Methylation_Status != "Unknown") %>%
  filter(Comment != "Fail") %>%
  dplyr::select(c("Sample_ID","Pyro_Methylation_Status")) %>% unique()

Training <- column_to_rownames(Training,"Sample_ID")


X1 <- na.omit(Metsum_ROI_long_ISLAND[rownames(Training),])

Y1 <- as.data.frame(Training[rownames(X1),])

rownames(Y1) <- rownames(X1)

colnames(Y1) <- "Methylation"

Y1$Methylation <- as.factor(Y1$Methylation)

Y1 <- Y1$Methylation

# Make PCA plot 
PCA.MGMT <- pca(X1, ncomp = 10, center = F, scale = F)

# Plot the eigenvalues of first 10 PCs
plot(PCA.MGMT)

# PCA plot of firt 2 principal components
plotIndiv(PCA.MGMT, group = Y1, legend = T, ind.names = F, title = "PCA, comp 1 & 2")

# Make initial PLS-DA (full data)

PLSDA.MGMT <- plsda(X1,Y1, ncomp = 10)

plotIndiv(PLSDA.MGMT, comp = 1:2, group = Y1, ind.names = F, ellipse = T, legend = T)

# perform M-fold cross-validation
perf.PLSDA.MGMT <- perf(PLSDA.MGMT, validation = "Mfold", folds = 6, nrepeat = 50, progressBar = T, auc = T)

# plot the perf results
plot(perf.PLSDA.MGMT, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")

perf.PLSDA.MGMT$choice.ncomp

# Tune number of variables to use in PLS-DA
list.keepX <- c(1:10,  seq(15, 50, 5))

tune.sPLSDA.MGMT <- tune.splsda(X1, Y1, ncomp = 2, # calculate for first 4 components
                                 validation = 'Mfold',
                                 folds = 6, nrepeat = 50, # use repeated cross-validation
                                 dist = 'max.dist', # use max.dist measure
                                 measure = "BER", # use balanced error rate of dist measure
                                 test.keepX = list.keepX,
                                 cpus = 2)

plot(tune.sPLSDA.MGMT, col = color.jet(2))

tune.sPLSDA.MGMT$choice.keepX

####### sPLS-DA based on tuning

sPLSDA.MGMT <- splsda(X1,Y1,ncomp = 2, keepX = c(25,1))

plotIndiv(sPLSDA.MGMT,  group = Y1, ind.names = F, ellipse = T, legend = T, title = "sPLS-DA, 25 features")

legend=list(legend = levels(Y1), # set of classes
            col = unique(color.mixo(Y1)), # set of colours
            title = "Methylation", # legend title
            cex = 0.7) # legend size

cim <- cim(sPLSDA.MGMT, row.sideColors = color.mixo(Y1),legend = legend)

plotLoadings(sPLSDA.MGMT)

plotVar(sPLSDA.MGMT)

#### feature stability

perf.sPLSDA.MGMT <- perf(sPLSDA.MGMT, 
                          folds = 6, nrepeat = 50, # use repeated cross-validation
                          validation = "Mfold", dist = "max.dist",  # use max.dist measure
                          progressBar = TRUE)

stable.comp1 <- perf.sPLSDA.MGMT$features$stable$comp1
barplot(stable.comp1, xlab = 'variables selected across CV folds', 
        ylab = 'Stability frequency',
        main = 'CpG feature stability')

#### testing sensitivity and accuracy with randomly splitting data into training and testing

train <- sample(1:nrow(X1), 100) # randomly select 50 samples in training
test <- setdiff(1:nrow(X1), train) # rest is part of the test set

X.train <- X1[train, ]
X.test <- X1[test,]
Y.train <- Y1[train]
Y.test <- Y1[test]

# train the model
train.MGMT <- splsda(X.train, Y.train, ncomp = 1, keepX = 25)

# predict the testset
predict.MGMT <- predict(train.MGMT, X.test, 
                                dist = "mahalanobis.dist")

# evaluate the prediction accuracy for the first two components
predict.comp1 <- predict.MGMT$class$mahalanobis.dist[,1]
table(factor(predict.comp1, levels = levels(Y1)), Y.test)


######## Predict the DenStem Samples

Testing <- MGMT_RunSum %>% 
  filter(Series == "DenStem") %>% 
  filter(Pyro_Methylation_Status != "Unknown") %>%
  filter(Comment != "Fail") %>%
  dplyr::select(c("Sample_ID","Pyro_Methylation_Status")) %>% unique()

Testing <- column_to_rownames(Testing,"Sample_ID")


X2 <- na.omit(Metsum_ROI_long_ISLAND[rownames(Testing),])

Y2 <- as.data.frame(Testing[rownames(X2),])

rownames(Y2) <- rownames(X2)

colnames(Y2) <- "Methylation"

Y2$Methylation <- as.factor(Y2$Methylation)

Y2 <- Y2$Methylation

#######
predict.MGMT <- predict(sPLSDA.MGMT, X2, 
                        dist = "mahalanobis.dist")
predict.comp2 <- predict.MGMT$class$mahalanobis.dist[,2]
table(factor(predict.comp2, levels = levels(Y2)), Y2)

Comparison <- as.data.frame(predict.comp2)
Comparison$truth <- Y2
