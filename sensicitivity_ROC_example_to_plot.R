### load data
load("sensicitivity_ROC_example_to_plot.RData")

### boxplot plot
library(ggplot2)
library(plotROC)  ### insatll.packages()


#### boxplot 1
ggplot(sensitivity_result2, aes(Method,data, fill= diseaseID))+geom_boxplot()

#### boxplot 2
ggplot(sensitivity_result2, aes(Method,data, fill= diseaseID))+ geom_boxplot()+ facet_wrap(~diseaseID)+ theme(axis.text.x=element_text(angle=90))

#### boxplot 3
ggplot(sensitivity_result2, aes(Method,data, fill= diseaseID))+ geom_boxplot()+ facet_grid(.~diseaseID)+ theme(axis.text.x=element_text(angle=90))

#### ROC
# compute AUC value
# roc_ssgsea <- roc(data$Type, data$ssgsea)
# print(AucHT_TG <- auc(roc_ssgsea))
# roc_gsva <- roc(data$Type, data$gsva)
# print(AucHT_TPO <- auc(rocHT_TPO))

# ROC-one line 
basicplot <- ggplot(results_SenSpe, aes(d = Type,m = Value, color = Method))+ geom_roc()
styleplot <- basicplot + style_roc(theme = theme_grey, xlab = "Specificity",ylab = "Sensitivity") + ggtitle("ROC")
styleplot

#### add annotation
p <- ggplot(results_SenSpe, aes(d = Type,m = Value, color = Method))+ 
  geom_roc() + 
  style_roc(theme = theme_grey, xlab = "Specificity",ylab = "Sensitivity") + ggtitle("ROC") +
  ggsci::scale_color_lancet()
auc <- calc_auc(p)  ## calculate the AUC
head(auc)  ## CHECK THE VALUE
p+annotate("text",x = .75, y = .45, 
           label = paste("AUC of ssgsea =", round(calc_auc(p)$AUC[3], 2))) +
  annotate("text",x = .75, y = .35, 
           label=paste("AUC of gsva =", round(calc_auc(p)$AUC[1], 2))) +
  annotate("text",x = .75, y = .25, 
           label = paste("AUC of plage =", round(calc_auc(p)$AUC[2], 2))) +
  annotate("text",x = .75, y = .15, 
           label=paste("AUC of zscore =", round(calc_auc(p)$AUC[4], 2)))

