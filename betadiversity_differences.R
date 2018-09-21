require(vegan)
require(ggplot2)

### Example dataset
# set.seed(123)
# t0 <- data.frame(x = rpois(n = 20,lambda = 5), y = rpois(n = 20,lambda = 5),timepoint = "t0",treat = rep(c("trattato","placebo"),each = 10))
# t1trattati <- data.frame(x = rpois(n = 10,lambda = 10), y = rpois(n = 10,lambda = 1),timepoint = "t1",treat = "trattato")
# t1placebo <- data.frame(x = rpois(n = 10,lambda = 5), y = rpois(n = 10,lambda = 5),timepoint = "t1",treat = "placebo")
# t0$sample <- paste("sample_",1:20,sep = "")
# t1trattati$sample <- paste("sample_",1:10,sep = "")
# t1placebo$sample <- paste("sample_",11:20,sep = "")
# 
# df <- rbind(t0,t1placebo,t1trattati)
# ggplot(df,aes(x,y,shape = timepoint,color = treat)) + geom_point() + geom_path(aes(group = sample))

distbeta <- function(beta, # Matrice della beta diversity
                     comp = c(1,2), # vettore di lunghezza 2 con le componenti su cui si vuole calcolare la distanza
                     metadata, # sample_data dell'oggetto su cui è stata calcolata la beta diversità
                     samplevar = NULL, # nome della variabile contenente la lista dei campioni
                     condvar, # nomi delli variabili su cui eseguire un test
                     timevar, # nome della variabile temporale
                     timediff = c("t0-t1")) # Istanti teporali su cui si è interessati a calcolare la distanza
{
  betadist <- beta # Matrice delle coordinate
  if(mean(rownames(betadist)==rownames(metadata))!=1)
    return(cat("I campioni non sono ordinati allo stesso modo nella matrice delle beta-diversity e nei metadati"))
  out <- NULL
  for(timepoint in timediff)
  {
    istanti <- unlist(strsplit(timepoint,"-"))
    index.tempo1 <- which(metadata[,timevar]==istanti[1])
    index.tempo2 <- which(metadata[,timevar]==istanti[2])
    betadist.tempo1 <- betadist[index.tempo1,as.vector(comp)]
    betadist.tempo2 <- betadist[index.tempo2,comp]
    
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

distbeta.beta <- function(beta, # Matrice della beta diversity
                     metadata, # sample_data dell'oggetto su cui è stata calcolata la beta diversità
                     samplevar = NULL, # nome della variabile contenente la lista dei campioni
                     condvar, # nomi delli variabili su cui eseguire un test
                     timevar, # nome della variabile temporale
                     timediff = c("t0-t1")) # Istanti teporali su cui si è interessati a calcolare la distanza
{
  betadist <- beta # Matrice delle coordinate
  if(mean(rownames(betadist)==rownames(metadata))!=1)
    return(cat("I campioni non sono ordinati allo stesso modo nella matrice delle beta-diversity e nei metadati"))
  out <- NULL
  for(timepoint in timediff)
  {
    istanti <- unlist(strsplit(timepoint,"-"))
    index.tempo1 <- which(metadata[,timevar]==istanti[1])
    index.tempo2 <- which(metadata[,timevar]==istanti[2])
    
    if(!is.null(samplevar)){
      coords <- apply(metadata[,samplevar],1,function(id) which(metadata[,samplevar] == id))
      dist <- lapply(coords,function(xy){
        if(length(xy)>1)
        {
          beta[xy[1],xy[2]]
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
        ord = match(campioni.t2,campioni.t1)
        dist(ord)
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

# distbeta(beta = df[,1:2],comp = 1:2,metadata = df[,3:5],condvar = "treat",samplevar = NULL,timevar = "timepoint",timediff = c("t0-t1"))

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

