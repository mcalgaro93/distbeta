require(vegan)
require(ggplot2)
require(phyloseq)

distbeta <- function(beta, # PCA matrix for beta diversity
                     comp = c(1,2), # Which PCA components do you want?
                     metadata, # as.data.frame(sample_data(phyloseq_object))
                     samplevar = NULL, # sample names variable
                     condvar, # interesting variable such as Treatment, Gender, BMI...
                     timevar, # time variable names
                     timediff = c("t0-t1")) # Timepoint differences
{
  betadist <- beta # coordinates
  if(mean(rownames(betadist)==rownames(metadata))!=1)
    return(cat("Different order of the samples between metadata and PCA coordinates"))
  out <- NULL
  for(timepoint in timediff)
  {
    istanti <- unlist(strsplit(timepoint,"-"))
    index.tempo1 <- which(metadata[,timevar]==istanti[1])
    index.tempo2 <- which(metadata[,timevar]==istanti[2])
    betadist.tempo1 <- betadist[index.tempo1,as.vector(comp)]
    betadist.tempo2 <- betadist[index.tempo2,as.vector(comp)]
    
    if(!is.null(samplevar)){
      campioni.t1 <- metadata[index.tempo1,samplevar]
      campioni.t2 <- metadata[index.tempo2,samplevar]
      
      if(length(campioni.t1)>length(campioni.t2)){
        betadist.tempo1 <- betadist.tempo1[match(campioni.t2,campioni.t1),]
      } else{
        betadist.tempo2 <- betadist.tempo2[match(campioni.t1,campioni.t2),]
      }
    }
    
    dist <- apply(betadist.tempo2-betadist.tempo1,1,function(x) sqrt(sum(x^2)))
    condition <- cbind(metadata[match(names(dist),rownames(metadata)),c(samplevar,condvar)],timepoint)
    dataframe <- data.frame(condition,dist = dist)
    out <- c(out,list(dataframe))
  }
  names(out) <- timediff
  class(out) <- "distbeta"
  return(out)
}

distbeta.beta <- function(beta, # beta diversity matrix
                     metadata, # as.data.frame(sample_data(phyloseq_object))
                     samplevar = NULL, # sample names variable
                     condvar, # interesting variable such as Treatment, Gender, BMI...
                     timevar, # time variable names
                     timediff = c("t0-t1")) # Timepoint differences
{
  betadist <- beta # betadiversity
  if(mean(rownames(betadist)==rownames(metadata))!=1)
    return(cat("Different order of the samples between metadata and betadiversity"))
  out <- NULL
  for(timepoint in timediff)
  {
    istanti <- unlist(strsplit(timepoint,"-"))
    index.tempo1 <- which(metadata[,timevar]==istanti[1])
    index.tempo2 <- which(metadata[,timevar]==istanti[2])
    
    if(!is.null(samplevar)){
      coords <- apply(metadata[,c(samplevar)],1,function(id) cbind(metadata[which(metadata[,samplevar] == id[1]),timevar],coord = which(metadata[,samplevar] == id[1])))
      
      dist <- lapply(coords,function(xy){
        if(sum(istanti %in% xy[,1])==2)
        {
          beta[xy[which(xy[,1] == istanti[1]),2],xy[which(xy[,1] == istanti[2]),2]]
        } else NULL
      })
      campioni.t1 <- as.matrix(metadata[index.tempo1,samplevar])
      campioni.t2 <- as.matrix(metadata[index.tempo2,samplevar])
      if(length(campioni.t1)>length(campioni.t2)){
        dist <- dist[index.tempo1]
        ord = match(campioni.t2,campioni.t1)
        dist <- dist[ord]
      } else{
        dist <- dist[index.tempo2]
        ord = match(campioni.t1,campioni.t2)
        dist <- dist[ord]
      }
    }
    
    condition <- cbind(metadata[match(names(dist),rownames(metadata)),c(samplevar,condvar)],timepoint)
    dataframe <- data.frame(condition,dist = unlist(dist))
    out <- c(out,list(dataframe))
  }
  names(out) <- timediff
  class(out) <- "distbeta"
  return(out)
}

#prova <- distbeta.beta(beta = beta,metadata = as.data.frame(sample_data(ps0)),samplevar = "ID",condvar = c("Treatment","Disease","Gender"),timevar = "TimePoint",timediff = "T0-T1")
#prova <- distbeta(beta = beta$vectors,comp = c(1,2),metadata = sampledf,condvar = c("Treatment","Disease","Gender","Age"),samplevar = "ID",timevar = "TimePoint",timediff = c("T0-T1"))         

ggplot.distbeta <- function(distbeta,timediff = NULL,var.formula = NULL,...) # Oggetto di tipo distbeta
{
  if(!is.null(var.formula) & !is.null(timediff)){
    sapply(distbeta[timediff], function(data){
      var = attr(terms.formula(var.formula), "term.labels")
      sapply(var, function(condition){
        if(length(unlist(strsplit(condition,":")))>1)
          data[,condition] <- apply(data[,unlist(strsplit(condition,":"))],1,paste,collapse = ":")
        print(ggplot(data,aes(x = data[,condition],y = dist, color = data[,condition])) + 
                geom_boxplot() + geom_jitter() + 
                ggtitle(paste(condition,levels(data$timepoint))) + 
                xlab(condition) + 
                labs(color = condition,y = "Distance") + 
                theme(axis.text.x = element_text(angle = 45,hjust = 1)))
      })
    })
  }
}

#prova.plot <- ggplot.distbeta(prova,timediff = c("T0-T1"),var = ~Treatment*Disease*Gender)

test.distbeta <- function(distbeta = NULL,var.formula,type = c(""))
{
  for(data in distbeta){
    data <- as.data.frame(data)
    if("" %in% type){
      m <- step(lm(as.formula(paste0("dist",paste0(var.formula,collapse = ""),collapse = "")),data = data),direction = "both")
      print(summary(m))
    }
  }
}

#test.distbeta(distbeta = prova,var.formula = ~ Gender*Disease*Treatment)

