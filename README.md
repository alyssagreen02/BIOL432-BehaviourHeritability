# BIOL432-BehaviourHeritability

## Team 5: Alyssa Green, Cameron Forsyth, Elliot Gavrin, Marco Yip, Eva Burguete-Innocente
[Link to dataset](https://doi.org/10.5061/dryad.b38k42m)

##Introduction

This project analyses meta-analysis data from "The Heritability of Behavior: A Meta-analysis" by Dochtermann et al. (2019). They evaluated the heritability of certain behaviours across taxa and based on several different traits. 

##Data

*Estimate.ID: unique ID for each estimate. 
*Study.ID: unique ID associated with the source paper. 
*Vert: vertebrate or invertebrate. 
*Thermy: endo- or ectotherm. 
*Phylum, Class, Order, Family, Genus, species.epithet: Phylogenetic information for the species associated with the estimate. 
*phylo: species name as it occurs in the phylogenetic tree. 
*D.W.: domesticated or wild species. 
*L.F.: estimate was made in lab or field conditions. 
*Behavior: behavioral category corresponding to estimate.
*N: sample size at the highest hierarchical level. 
*N.C.: alternative sample size to address whether inflated precision biased results. 
*Heritability: h2 estimate. 
*Zr: transformed estimate. Transformed using Fisher's z-transformation by the original authors. 

##Analyses

1. Downloaded the csv file from Dochtermann et al. 
2. Produced phylogenetic tree of all taxa included in the dataset.
3. Ran linear mixed effects models to evaluate the effects of taxonomic rank and traits and heritability. 
4. Ran linear mixed effects models to evaluate the effects of different behaviours on heritability. 
