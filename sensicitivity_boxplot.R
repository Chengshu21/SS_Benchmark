### load data
load("sensicitivity_example_to_plot.RData")

### boxplot plot

library(ggplot2)

####boxplot 1
ggplot(sensitivity_result2, aes(Method,data, fill= diseaseID))+geom_boxplot()

####boxplot 2
ggplot(sensitivity_result2, aes(Method,data, fill= diseaseID))+ geom_boxplot()+ facet_wrap(~diseaseID)+ theme(axis.text.x=element_text(angle=90))

####boxplot 3
ggplot(sensitivity_result2, aes(Method,data, fill= diseaseID))+ geom_boxplot()+ facet_grid(.~diseaseID)+ theme(axis.text.x=element_text(angle=90))

