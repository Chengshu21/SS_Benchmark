#### select NON-significant pathways of each cluster at first
##################################
#### get TN and FN

get_tn.fn = function(nonsig.pathway_result, target_pathway){
  tn_count <- 0
  fn_count <- 0
  for (i in 1:length(rownames(nonsig.pathway_result))) 
  {
    for (j in 1:length(target_pathway))
    {
      if (rownames(nonsig.pathway_result)[i] == as.character(target_pathway[j])) {
        fn_count <- fn_count + 1
        FN[[fn_count]] <- target_pathway[j]
      }
      if (rownames(nonsig.pathway_result)[i] != as.character(target_pathway[j])) {
        tn_count <- tn_count + 1
        TN[[tn_count]]<- rownames(nonsig.pathway_result)[i]
      }
    }
  }
  result = list("tn_count" = nrow(nonsig.pathway_result)-fn_count,
                "fn_count" = fn_count,
                "FN" = FN)
}

