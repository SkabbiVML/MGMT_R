df <- FOUR_CpGs_Group %>% select(Average_Nano,Known_status)

df$Known_status <- ifelse(df$Known_status=="Methylated",1,0)

glm.fit <- glm(df$Known_status ~ df$Average_Nano, family = binomial)

library(pROC)

par(pty = "s")


roc(df$Known_status, 
    glm.fit$fitted.values, 
    legacy.axes = T, 
    plot = T, 
    percent = T, 
    xlab ="False Positive percentage",
    ylab = "True Positive percentage",
    lwd = 3,
    print.auc = T)

roc.info<- roc(df$Known_status, glm.fit$fitted.values, legacy.axes = T)
roc.df <- data.frame(tpp=roc.info$sensitivities*100, 
                     fpp=(1 - roc.info$specificities)*100,
                     thresholds=roc.info$thresholds)

par(pty = "m")

ggplot(roc.df)