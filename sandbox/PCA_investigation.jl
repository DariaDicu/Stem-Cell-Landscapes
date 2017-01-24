#This script aims to consider the viability of using Principal
#Component Analysis (PCA) as a dimensionality reduction
#technique to aid in visualising our gene expression data.

#Supporting documentation found at:
#http://multivariatestatsjl.readthedocs.io/en/latest/pca.html

using MultivariateStats

# 'One can use the fit method to perform PCA over a given dataset'
# Perform PCA over the data given in a matrix, X, where each column
# of X is an observation. fit(PCA, X; ...)'

#Let (d,n)=size(X) be respectively the input dimension and the
#number of observations

X = rand(10, 30)

PCA_out = fit(PCA, X; maxoutdim = 2)

plot(PCA_out)

comps = transform(PCA_out, X)

plot(plot(comps[:,1]), plot(comps[:,2]), plot(comps[:,1], comps[:,2]))
