setwd("~/MGMT_stuff") # Home desktop


colnames(Metsum_ROI_long) <- c(-7:-1, 1:98,c("+1","+2", "+3","+4", "+5", "+6", "+7", "+8", "+9", "+10"))

X1 <- Metsum_ROI_long[c(2:129),]

#X1 <- logratio.transfo(X1,logratio = "CLR")

Y1 <- MGMT_RunSum %>% 
  filter(Series != "DenStem") %>% 
  dplyr::select(SampleID,Known_status) %>% 
  distinct() %>% 
  na.omit()

rownames(Y1) <- Y1$SampleID 

Y1 <- Y1[rownames(X1),2]


library(mixOmics)

#### SImple PCA

pca.meth <- pca(X1, ncomp = 10,center = TRUE, scale = TRUE)

plotIndiv(pca.meth, 
          ind.names = F,
          group = Y1 )

spca.meth <- spca(X1, ncomp = 10, center = TRUE)

plotIndiv(spca.meth, 
          ind.names = F,
          group = Y1 )

plot(pca.meth)



##### Simple PLS-DA

plsda.Meth <- plsda(X1,Y1, ncomp = 5)

plotIndiv(plsda.Meth, 
          ind.names = F,
          group = Y1 )


######### sPLS-DA


splsda.Meth <- splsda(X1,Y1,ncomp = 5)


plotIndiv(splsda.Meth, comp = 1:2,
          group = Y1,
          ind.names = F,
          ellipse = T)



perf.plsda.Meth <- perf(plsda.Meth, validation = 'Mfold', folds = 5, 
                          progressBar = T,  # Set to TRUE to track progress
                          nrepeat = 50)         # We suggest nrepeat = 50

plot(perf.plsda.Meth)

perf.plsda.Meth$choice.ncomp# Select 2 components, max distance

# grid of possible keepX values that will be tested for each component
list.keepX <- c(1:10,  seq(12, 30, 2))

# undergo the tuning process to determine the optimal number of variables
tune.splsda.Meth <- tune.splsda(X1, Y1, ncomp = 2, # calculate for first 4 components
                                 validation = 'Mfold',
                                 folds = 5, nrepeat = 50, # use repeated cross-validation
                                 dist = 'max.dist', # use max.dist measure
                                 measure = "BER", # use balanced error rate of dist measure
                                 test.keepX = list.keepX,
                                 progressBar = T,
                                 cpus = 4) # allow for paralleliation to decrease runtime


plot(tune.splsda.Meth)


tune.splsda.Meth$choice.ncomp$ncomp # 2 components

tune.splsda.Meth$choice.keepX # comp 1: 30, comp2:10

# store values for final model

optimal.ncomp <- tune.splsda.Meth$choice.ncomp$ncomp
optimal.keepX <- tune.splsda.Meth$choice.keepX[1:optimal.ncomp]

#### Final spls-da model

final.splsda.Meth <- splsda(X1,Y1,ncomp = optimal.ncomp, scale = T, 
                            keepX = optimal.keepX)


#### Final spls-Da plot

plotIndiv(final.splsda.Meth, comp = c(1,2), group = Y1, ind.names = F, ellipse = T, ellipse.level = 0.90)


######## heatmap

legend=list(legend = as.factor(unique(MGMT_RunSum$Known_status)), # 
            col = RColorBrewer::brewer.pal(3,"Set1"), # set of colours
            title = "Known methylation", # legend title
            cex = 0.7) # legend size


cim(final.splsda.Meth, legend = legend, row.sideColors = RColorBrewer::brewer.pal(3,"Set1"))


final.plsda.Rad <- plsda(X1,Y1, ncomp = 2) 


plotIndiv(final.plsda.Rad)

perf.final.plsda.Rad <- perf(final.plsda.Rad, validation = 'Mfold', 
                                folds = 5, 
                                progressBar = T, # TRUE to track progress
                                nrepeat = 100) 

plot(perf.final.plsda.Rad, sd = TRUE, legend.position = 'horizontal') #teo components, max.dist

perf.final.plsda.Rad$error.rate$BER[, 'max.dist']

perf.final.plsda.Rad$error.rate$overall[, 'max.dist']

perf.final.plsda.Rad$error.rate.class$max.dist

background.max.Rad <- background.predict(final.plsda.Rad, 
                                            comp.predicted = 2,
                                            dist = 'max.dist')

plotIndiv(final.plsda.Rad, comp = 1:2,
          ind.names = FALSE, title = 'PLS-DA Pyro evaluated methylation',
          legend = TRUE,  background = background.max.Rad)

list.keepX <- c(seq(5,100,5))

tune.splsda.Rad <- tune.splsda(X1, Y1, ncomp = 2,
                                  validation = 'Mfold', 
                                  folds = 5 , dist = 'max.dist', 
                                  test.keepX = list.keepX, nrepeat = 100,
                                  progressBar = T)

head(tune.splsda.Rad$error.rate)

plot(tune.splsda.Rad, sd = TRUE)

select.keepX <- tune.splsda.Rad$choice.keepX[1:2]

# sparse PLS-DA models
splsda.Rad <- splsda(X1, Y1, ncomp = 2, keepX = 10)

perf.splsda.Rad <- perf(splsda.Rad, folds = 5, validation = "Mfold", 
                           dist = "max.dist", progressBar = T, nrepeat = 100)

plotIndiv(splsda.Rad, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - Radium samples", legend.title = 'PyroSeq')

auc.plsda <- auroc(splsda.Rad)

plotLoadings(splsda.Rad, contrib = 'max', method = 'man', title = "Loadings comp 1", legend = F, ndisplay = 40)

legend=list(legend = levels(Y1), # set of classes
            col = unique(color.mixo(Y1)), # set of colours
            title = "Methylation", # legend title
            cex = 0.7) # legend size

cim <- cim(splsda.Rad, 
           legend = legend)

# use the model on the Xtest set

VML <- read.csv("MGMT_VML/VML_MGMT_Methylation_summary.csv", header = T, row.names = 1)
colnames(VML) <- c(-6:-1, 1:98,c("+1","+2", "+3","+4", "+5", "+6", "+7", "+8", "+9"))

X2 <- VML[c(6:nrow(VML)),]

VML_names <- rownames(X2)

X2 <- sapply(X2, as.numeric)
rownames(X2) <- VML_names

Anno_VML <- row_anno[rownames(X2),]

library(dplyr)
Y2 <-Anno_VML$PyroSeq
#Y2$PyroSeq <- as.factor(Y2$PyroSeq)

# Predict methylation status of VML samples based on model built on Radium samples
predict.splsda.VML <- predict(splsda.Rad, X2, 
                                dist = "max.dist")
predict.comp2 <- predict.splsda.VML$class$max.dist[,2]
table(factor(predict.comp2, levels = unique(row_anno$PyroSeq)), Y2)

Comparison <- cbind(predict.comp2,Y2)

####### PLS-DA VML

Y2 <- Anno_VML$PyroSeq

plsda.VML <- plsda(X2,Y2, ncomp = 2)

plotIndiv(plsda.VML, ind.names = T, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - VML samples", legend.title = 'PyroSeq')

######## PCA

tune.pca.Rad <- tune.pca(X1, ncomp = 5, scale = TRUE)
plot(tune.pca.Rad, main = "Screeplot of PCA performance, Radium")
tune.pca.Rad$cum.var       # Outputs cumulative proportion of variance
plot(tune.pca.Rad$cum.var)

## ----echo=TRUE, message=FALSE----------------------------------------------------------------------------------------------------------------------
final.pca.Rad <- pca(X1, ncomp = 3, center = TRUE, scale = TRUE)

plotIndiv(final.pca.Rad,
          ind.names = F, # Show row names of samples
          group = Anno_Radium$PyroSeq,
          
          title = 'PCA Radium samples',
          legend = TRUE, legend.title = 'PyroSeq')

plotLoadings(final.pca.Rad)

final.pca.VML <- pca(X2, ncomp = 3, center = TRUE, scale = TRUE)

plotIndiv(final.pca.VML,
          ind.names = T, # Show row names of samples
          group = Anno_VML$NPexon1,
          title = 'Pyroseq status VML, PCA comp 1 - 2',
          legend = TRUE, legend.title = 'PyroSeq')

