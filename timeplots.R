library(phyloseq)
library(vegan)
library(ggplot2)

dist_matrix <- phyloseq::distance(phyloseq_object , method = "unifrac")
beta <- as.matrix(dist_matrix)
db <- distbeta.beta(beta = beta,metadata = as.data.frame(sample_data(phyloseq_object)),samplevar = "Sample",condvar = c("Treatment","Gender","Age"),timevar = "TimePoint",timediff = c("T0-T1","T1-T2","T2-T3"))
df <- rbind(db$`T0-T1`,db$`T1-T2`,db$`T2-T3`)

ggplot(df,aes(x = timepoint,y = dist,color = Sample)) + 
  geom_point() + 
  geom_path(aes(group = Sample),size = 1) +
  facet_grid(Gender~Treatment) 
