topFeaturesOmnibus <-
function (models, contrasts, adjust.method = "BH", sort = TRUE,
    alpha = 1)
{

    qrL <- qr(contrasts)
    Lh <- qr.Q(qrL)[,1:qrL$rank]
    logFCs <- apply(contrasts, 2, function(models, L) sapply(models, getContrast, L = L), models = models)
    Fstat <- sapply(models,
        function(model, Lh, r){
          if (getFitMethod(model)=="fitError")
            return(NA) else {
              est <- getCoef(model)%*%Lh
              varEst <- t(Lh)%*%(getVcovUnscaled(model) * getVarPosterior(model))%*%Lh
              Fstat <- est %*% solve(varEst) %*% t(est) / r
              return(Fstat)
            }
          },
          Lh = Lh,
          r = qrL$rank)
    df <- sapply(models, getDfPosterior)
    pval <- pf(Fstat, qrL$rank, df, lower.tail=FALSE)
    adjPval <- p.adjust(pval, method = adjust.method)
    out <- data.frame(logFCs, df1 = qrL$rank, df2 = df, Fstat, pval, adjPval)
    if (sort) {
        if (alpha < 1) {
            ids <- adjPval < alpha
            out <- na.exclude(out[ids, ])
        }
        return(out[order(out$pval), ])
    }
    else {
        return(out)
    }
}

setGeneric("omnibusTest",function(object,...) standardGeneric("omnibusTest"))

setMethod("omnibusTest","QFeatures",
          function(object,
                   i,
                   contrasts,
                   adjust.method="BH",
                   modelColumn="msqrobModels",
                   resultsColumnName="omnibusTest",
                   overwrite=FALSE){
            if (is.null(object[[i]])) stop(paste0("QFeatures object does not contain an assay with the name ",i))
            if(!(modelColumn %in% colnames(rowData(object[[i]])))) stop(paste0("There is no column named \'", modelColumn,"\' with stored models of an msqrob fit in the rowData of assay ",i,"of the QFeatures object."))
            if(is.null(colnames(contrasts)) & resultsColumnName=="") resultsColumnName<-"omnibusTest"
            if((resultsColumnName %in% colnames(rowData(object[[i]])))&!overwrite) stop(paste0("There is/are already column(s) named \'", resultsColumnName, "\' in the rowData of assay ",i," of the QFeatures object, set the argument overwrite=TRUE to replace the column(s) with the new results or use another name for the argument resultsColumnName"))
            rowData(object[[i]])[[resultsColumnName]] <- topFeaturesOmnibus(rowData(object[[i]])[,modelColumn], contrasts = contrasts, adjust.method = adjust.method, sort = FALSE, alpha = 1)
            return(object)
            })
