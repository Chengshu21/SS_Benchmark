#### select significant pathways of each cluster at first
##################################
#### get TP and FP
TP <- list()
TN <- list()
FP <- list()
FN <- list()

get_tp.fp = function(sig.pathway_result, target_pathway){
  tp_count <- 0
  fp_count <- 0
  TP <- list()
  for (i in 1:length(rownames(sig.pathway_result))) 
  {
    for (j in 1:length(target_pathway))
    {
      if (rownames(sig.pathway_result)[i] == as.character(target_pathway[j])) {
        tp_count <- tp_count + 1
        TP[[tp_count]]<- target_pathway[j]
      }
      if (rownames(sig.pathway_result)[i] != as.character(target_pathway[j])) {
        fp_count <- fp_count + 1
        FP[[fp_count]] <- rownames(sig.pathway_result)[i]
      }
    }
  }
  result = list("tp_count" = tp_count,
                "TP" = TP, 
                "fp_count" = nrow(sig.pathway_result)-tp_count)
}

