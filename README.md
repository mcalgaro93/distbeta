# distbeta
Extract distances for each sample between different timepoint starting from a dissimilarity matrix.

Imagine a dissimilarity matrix containing beta diversities for a follow-up study. Several time points exist and we want to study beta diversity differences for each individual, between couple of timepoints. 

If treatment is responsible for a change in host microbiome, we'll maybe observe bigger changes on treated individuals beta-diversity between subsequent time points. 

There are several functions:
- distbeta computes differences between samples using PCA coordinates, the user can choose how many components to use. Metadata rownames must be in the same order of input matrix. The user has to specify which metadata column names contain sample names, useful variable and timepoint variable. Timepoints differences are computed for specified timepoints only (eg. "T0-T1","T1-T2"...);
- distbeta.beta computes differences between samples using the dissimilarity matrix. It is more complete than distbeta. There is no need of PCA components, the input matrix is the entire #sample x #sample dissimilarity matrix;
- ggplot.distbeta plots as many graphics as defined in timediff and var.formula;
- test.distbeta, in a very primitive way, searches the most parsimonious model through stepwise selection beginning from var.formula. Significant variables are those who might be associated with distances.

# timeplots
In order to obtain time-series plot I suggest the creation of a data.frame object like explained in timeplots.R
