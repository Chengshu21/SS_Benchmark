#######gene sets(pathways) file######
# geneset file
### 1.genesetcollection object used in "gsva" package
library(GSVA)

KEGG_genesetcollection = getGmt("F:/lab_data/3. MSigDB files/c2.all.v6.2.entrez.gmt",
                                geneIdType = EntrezIdentifier(),
                                collectionType = BroadCollection(category="c2"))
KEGG_genesetcollection = KEGG_genesetcollection[grep("^KEGG", names(KEGG_genesetcollection))]
# grep("^REACTOME", names(genesetcollection)),
# grep("^BIOCARTA", names(genesetcollection)))]
KEGG_genesetcollection

### 2.geneset list object
read.gmt = function(file){
  if(!grepl("\\.gmt$",file)[1]){stop("Pathway information must be a .gmt file")}
  geneSetDB = readLines(file)                                ##read in the gmt file as a vector of lines
  geneSetDB = strsplit(geneSetDB,"\t")                       ##convert from vector of strings to a list
  names(geneSetDB) = sapply(geneSetDB,"[",1)                 ##move the names column as the names of the list
  geneSetDB = lapply(geneSetDB, "[",-1:-2)                   ##remove name and description columns
  geneSetDB = lapply(geneSetDB, function(x){x[which(x!="")]})##remove empty strings
  return(geneSetDB)
}
KEGG_genesetlist = read.gmt("F:/lab_data/3. MSigDB files/c2.all.v6.2.entrez.gmt")
KEGG_genesetlist = KEGG_genesetlist[grep("^KEGG", names(KEGG_genesetlist))]
# grep("^REACTOME", names(genesetlist)),
# grep("^BIOCARTA", names(genesetlist)))]
# KEGG_pathwaynameslist = as.matrix(KEGG_genesetlist)
# KEGG_pathwaynameslist = list(gs = KEGG_genesetlist, pathwaynames = as.list(names(KEGG_pathwaynameslist)))
# str(KEGG_genesetlist)

### 3. Pathway list used in "pathifier" package
KEGG_pathwaynamelist = as.character(do.call(rbind, as.list(names(KEGG_genesetlist))))
KEGG_pathwaynamelist

# Load Genesets annotation
#  Generate a list that contains genes in genesets
genesets = KEGG_genesetlist

# Generate a list that contains the names of the genesets used
H00014 = list("KEGG_NON_SMALL_LUNG_CANCER_H00014" = "KEGG_NON_SMALL_LUNG_CANCER",
              "KEGG_MAPK_SIGNALING_PATHWAY_H00014" = "KEGG_MAPK_SIGNALING_PATHWAY",
              "KEGG_ERBB_SIGNALING_PATHWAY_H00014" = "KEGG_ERBB_SIGNALING_PATHWAY",
              "KEGG_RAS_SIGNALING_PATHWAY_H00014" = "KEGG_RAS_SIGNALING_PATHWAY",
              "KEGG_CALCIUM_SIGNALING_PATHWAY_H00014" = "KEGG_CALCIUM_SIGNALING_PATHWAY",
              "KEGG_CELL_CYCLE_H00014" = "KEGG_CELL_CYCLE",
              "KEGG_P53_SIGNALING_PATHWAY_H00014" = "KEGG_P53_SIGNALING_PATHWAY",
              "KEGG_PI3K_AKT_SIGNALING_PATHWAY_H00014" = "KEGG_PI3K_AKT_SIGNALING_PATHWAY")
H01103 = list("KEGG_COMPLEMENT_AND_COAGULATION_CASCADES_H01103" = "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",
              "KEGG_INFLAMMATORY_MEDIATOR_REGULATION_OF_TRP_CHANNELS_H01103" = "KEGG_INFLAMMATORY_MEDIATOR_REGULATION_OF_TRP_CHANNELS",
              "KEGG_PLATELET_ACTIVATION_H01103" = "KEGG_PLATELET_ACTIVATION",
              "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY_H01103" = "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY")
H00079 = list("KEGG_ASTHMA_H00079" = "KEGG_ASTHMA",
              "KEGG_FC_EPSILON_RISIGNALING_PATHWAY_H00079" = "KEGG_FC_EPSILON_RISIGNALING_PATHWAY",
              "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY_H00079" = "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
              "KEGG_CELL_ADHESION_MOLECULES(CAMs)_H00079" = "KEGG_CELL_ADHESION_MOLECULES(CAMs)" ,
              "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY_H00079" = "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
              "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION_H00079" = "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
              "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION_H00079" = "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
              "KEGG_JAK_STAT_SIGNALING_PATHWAY_H00079" = "KEGG_JAK_STAT_SIGNALING_PATHWAY")
H00342 = list("KEGG_TUBERCULOSIS_H00342" = "KEGG_TUBERCULOSIS",
              "KEGG_APOPTOSIS_H00342" = "KEGG_APOPTOSIS",
              "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION_H00342" = "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
              "KEGG_JAK_STAT_SIGNALING_PATHWAY_H00342" = "KEGG_JAK_STAT_SIGNALING_PATHWAY",
              "KEGG_MAPK_SIGNALING_PATHWAY_H00342" = "KEGG_MAPK_SIGNALING_PATHWAY",
              "KEGG_TOLL_LIKERECEPTOR_SIGNALING_PATHWAY_H00342" = "KEGG_TOLL_LIKERECEPTOR_SIGNALING_PATHWAY",
              "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY_H00342" = "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY")
H01714 = list("KEGG_CITRATE_CYCLE_TCA_CYCLE_H01714" = "KEGG_CITRATE_CYCLE_TCA_CYCLE", 
              "KEGG_HEDGEHOG_SIGNALING_PATHWAY_H01714" = "KEGG_HEDGEHOG_SIGNALING_PATHWAY",
              "KEGG_WNT_SIGNALING_PATHWAY_H01714 " = "KEGG_WNT_SIGNALING_PATHWAY",
              "KEGG_PI3K_AKT_SIGNALING_PATHWAY_H01714" = "KEGG_PI3K_AKT_SIGNALING_PATHWAY", 
              "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY_H01714" = "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
              "KEGG_TGF_BETA_SIGNALING_PATHWAY_H01714" = "KEGG_TGF_BETA_SIGNALING_PATHWAY", 
              "KEGG_ECM_RECEPTOR_INTERACTION_H01714" = "KEGG_ECM_RECEPTOR_INTERACTION ",
              "KEGG_FOCAL_ADHESION_H01714" = "KEGG_FOCAL_ADHESION ",
              "KEGG_VEGF_SIGNALING_PATHWAY_H01714" = "KEGG_VEGF_SIGNALING_PATHWAY",
              "KEGG_SPHINGOLIPID_METABOLISM_H01714 " = "KEGG_SPHINGOLIPID_METABOLISM ",
              "KEGG_NOTCH_SIGNALING_PATHWAY_H01714" = "KEGG_NOTCH_SIGNALING_PATHWAY",
              "KEGG_ERBB_SIGNALING_PATHWAY_H01714" = "KEGG_ERBB_SIGNALING_PATHWAY",
              "KEGG_ABC_TRANSPORTERS_H01714" = "KEGG_ABC_TRANSPORTERS", 
              "KEGG_CELL_ADHESION_MOLECULES_CAMS_H01714" = "KEGG_CELL_ADHESION_MOLECULES_CAMS",
              "KEGG_CELL_CYCLE_H01714" = "KEGG_CELL_CYCLE ",
              "KEGG_OOCYTE_MEIOSIS_H01714 " = "KEGG_OOCYTE_MEIOSIS ",
              "KEGG_APOPTOSIS_H01714"= "KEGG_APOPTOSIS", 
              "KEGG_P53_SIGNALING_PATHWAY_H01714" = "KEGG_P53_SIGNALING_PATHWAY",
              "KEGG_HOMOLOGOUS_RECOMBINATION_H01714"= "KEGG_HOMOLOGOUS_RECOMBINATION",
              "KEGG_BASE_EXCISION_REPAIR_H01714" = "KEGG_BASE_EXCISION_REPAIR",
              "KEGG_NUCLEOTIDE_EXCISION_REPAIR_H01714" = "KEGG_NUCLEOTIDE_EXCISION_REPAIR",
              " KEGG_MISMATCH_REPAIR_H01714 " =" KEGG_MISMATCH_REPAIR", 
              "KEGG_NON_HOMOLOGOUS_END_JOINING_H01714"  = "KEGG_NON_HOMOLOGOUS_END_JOINING", 
              "KEGG_ENDOCYTOSIS_H01714" = "KEGG_ENDOCYTOSIS",
              "KEGG_LYSOSOME_H01714" = "KEGG_LYSOSOME",
              "KEGG_PEROXISOME_H01714 "="KEGG_PEROXISOME ",
              "KEGG_SNARE_INTERACTIONS_IN_VESICULAR_TRANSPORT_H01714" = "KEGG_SNARE_INTERACTIONS_IN_VESICULAR_TRANSPORT",
              "KEGG_VASCULAR_SMOOTH_MUSCLE_CONTRACTION_H01714" = "KEGG_VASCULAR_SMOOTH_MUSCLE_CONTRACTION", 
              "KEGG_REGULATION_OF_ACTIN_CYTOSKELETON_H01714 " = "KEGG_REGULATION_OF_ACTIN_CYTOSKELETON ",
              "KEGG_OXIDATIVE_PHOSPHORYLATION_H01714" = "KEGG_OXIDATIVE_PHOSPHORYLATION" ,
              "KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS_H01714"= "KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS",
              "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION_H01714" = "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
              "KEGG_ADIPOCYTOKINE_SIGNALING_PATHWAY_H01714" = "KEGG_ADIPOCYTOKINE_SIGNALING_PATHWAY",
              "KEGG_MAPK_SIGNALING_PATHWAY_H01714" = "KEGG_MAPK_SIGNALING_PATHWAY",
              "KEGG_TUBERCULOSIS_H01714" = "KEGG_TUBERCULOSIS",
              "KEGG_PENTOSE_PHOSPHATE_PATHWAY_H01714" = "KEGG_PENTOSE_PHOSPHATE_PATHWAY",
              "KEGG_HEMATOPOIETIC_CELL_LINEAGE_H01714"  = "KEGG_HEMATOPOIETIC_CELL_LINEAGE",
              "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES_H01714 " = "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES ",
              "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY_H01714" = "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
              "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY_H01714" = "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY" ,
              "KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY_H01714" = "KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY",
              "KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY_H01714"= "KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY",
              "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY_H01714" = "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",
              "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION_H01714" = "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
              "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY_H01714" = "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY ",
              "KEGG_FC_EPSILON_RI_SIGNALING_PATHWAY_H01714" = "KEGG_FC_EPSILON_RI_SIGNALING_PATHWAY",
              "KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS_H01714 " = "KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS",
              "KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION_H01714" = "KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION", 
              "KEGG_INTESTINAL_IMMUNE_NETWORK_FOR_IGA_PRODUCTION_H01714" = "KEGG_INTESTINAL_IMMUNE_NETWORK_FOR_IGA_PRODUCTION",
              "KEGG_CHEMOKINE_SIGNALING_PATHWAY_H01714" = "KEGG_CHEMOKINE_SIGNALING_PATHWAY")
Disease_pathwaynames = list(H00014,H01103,H00079,H00342,H01714)

####target.pathway

target.pathway_H00014 = Disease_pathwaynames[[1]]
target.pathway_H01103 = Disease_pathwaynames[[2]]
target.pathway_H00079 = Disease_pathwaynames[[3]]
target.pathway_H00342 = Disease_pathwaynames[[4]]
target.pathway_H01714 = Disease_pathwaynames[[5]]


str(target.pathway_H00014) ### H00342


#########################################
#####pathway data used in "pathifier"
pathwaynames = Disease_pathwaynames
# Generate a list that contains the previos two lists: genesets and their names
PATHWAYS = list()
PATHWAYS$gs = genesets
PATHWAYS$pathwaynames = pathwaynames

##############data preparation
normals_GSE10245 = as.vector(as.logical(GSE10245_pdata$Normal))
exp.matrix = as.matrix(GSE10245_setentrez)

####function to get new data to be used in"pathifier"

get_new_data = function(eSet, setentrez){
  
  # Prepare data and parameters 
  # Extract information from binary phenotypes. 1 = Control, 0 = Disease
  normals = as.vector(as.logical(pData(eSet)$Normal))
  exp_matrix = as.matrix(setentrez)
  
  # Calculate MIN_STD
  New_exp_matrix = exp_matrix[,as.logical(normals)]
  nsd = apply(New_exp_matrix, 1, sd)
  min_std = quantile(nsd, 0.25)
  
  # Calculate MIN_EXP
  min_exp = quantile(as.vector(exp_matrix), 0.1) # Percentile 10 of data
  
  # Filter low value genes. At least 10% of samples with values over min_exp
  # Set expression levels < MIN_EXP to MIN_EXP
  greater = apply(exp_matrix, 1, function(x) x > min_exp)
  new_greater = apply(greater, 2, mean)
  new_greater = names(new_greater)[new_greater > 0.1]
  exp_matrix = exp_matrix[new_greater,]
  exp_matrix[exp_matrix < min_exp] = min_exp
  
  # Set maximum 5000 genes with more variance
  V = names(sort(apply(exp_matrix, 1, var), decreasing = T))#[1:10000]
  V = V[!is.na(V)]
  exp_matrix = exp_matrix[V,]
  genes = rownames(exp_matrix) # Checking genes
  allgenes = as.vector(rownames(exp_matrix))

  # Generate a list that contains previous data: gene expression, normal status,
  # and name of genes
  DATASET = list()
  DATASET$allgenes = allgenes
  DATASET$normals = normals
  DATASET$data = exp_matrix
  DATASET$min_std = min_std
  DATASET$min_exp = min_exp
  
  result = DATASET
}

new_GSE10245data = get_new_data(GSE10245eSet, GSE10245_setentrez)
new_GSE106986data = get_new_data(GSE106986eSet, GSE106986_setentrez)
new_GSE1122data = get_new_data( GSE1122eSet,  GSE1122_setentrez)
new_GSE11906data = get_new_data(GSE11906eSet, GSE11906_setentrez)
new_GSE12472.1data = get_new_data(GSE12472.1eSet, GSE12472.1_setentrez)
new_GSE12472.2data = get_new_data(GSE12472.2eSet, GSE12472.2_setentrez)
new_GSE18842data = get_new_data(GSE18842eSet, GSE18842_setentrez)
new_GSE35571data = get_new_data(GSE35571eSet, GSE35571_setentrez)
new_GSE37768data = get_new_data(GSE37768eSet, GSE37768_setentrez)
new_GSE42057data = get_new_data(GSE42057eSet, GSE42057_setentrez)
new_GSE50834data = get_new_data(GSE50834eSet, GSE50834_setentrez)
new_GSE52819data = get_new_data(GSE52819eSet, GSE52819_setentrez)







