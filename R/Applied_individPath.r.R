applied_individPath = function(tumorFile, refFile, PathwayFile, cutoff){
    source(SRGgenePair.R)   ## check in "individPath" package
    source(individPathCal.R)  ## check in "individPath" package
    ControlData <- as.matrix(read.table(refFile, header=T, sep="\t"));
    CaseData <- as.matrix(read.table(tumorFile, header=T, sep="\t"));
    
    if(!is.matrix(ControlData)|!is.matrix(CaseData)){
        stop("tumorFile or refFile must be matrix!\n");
    }
    SampleInfo <- colnames(CaseData)[-1];
    NumSample <- length(SampleInfo);
    
    if (nrow(CaseData)!= nrow(ControlData)){
        stop("tumorFile and refFile must have the same number of rows!\n");
    }
    
    if (is.na(PathwayFile)){
        stop("Please input pathway data!\n");
    }
    PathData <- read.table(PathwayFile, header=F,sep="\t",fill=T);
    PathName <- as.matrix(PathData[,1]);
      
    ###---------Identifying Stable gene pair -------###
    print("Identifying stable and reversal intra-pathway gene pairs");
    GP.result <- SRGgenePair(ControlData, CaseData, PathData, cutoff);
    StableGP <- GP.result$BG.GenePairs;
    NumStable <- nrow(StableGP);
    ReversalGP <- GP.result$ReversalStat;
    PathGP <- GP.result$PathGP;
    
    ###---------Individualized altered pathway -------###
    Result <- NULL;
    AD.result <- NULL;
    for( i in 1:NumSample){
        print(paste("individPath_processing ", i, "/", NumSample," : ", SampleInfo[i], sep=""));
        patient <- CaseData[,i+1];
        names(patient) <- CaseData[,1];
        
        Result.tmp <- individPathCal (patient, StableGP, ReversalGP[,i], NumStable, PathGP);        
        Result <- cbind(Result,Result.tmp);
        AD.p <- as.matrix(p.adjust(Result.tmp, "BH"));
        AD.result <- cbind(AD.result,AD.p);
    }
    Last.Result <- cbind( PathName, Result);
    Last.ADResult <- cbind( PathName, AD.result);
    
    colnames(Last.Result)[2:(NumSample+1)] <- SampleInfo;
    colnames(Last.Result)[1] <- "PathwayID";
    colnames(Last.ADResult)[2:(NumSample+1)] <- SampleInfo;
    colnames(Last.ADResult)[1] <- "PathwayID";
    
    result = list("BH_Result" = Last.Result,"pvalue_Result" = Last.ADResult)
    result
}
