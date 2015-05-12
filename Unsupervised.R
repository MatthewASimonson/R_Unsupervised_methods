
# Unsupervised Methods:

########################
# CLUSTERING EXAMPLES:
########################

setwd("~/Desktop/Data_Science/zmPDSwR-master/Protein")
protein <- read.table(“protein.txt”, sep=”\t”, header=TRUE)
summary(protein)

# Title: Rescaling the dataset

vars.to.use <- colnames(protein)[-1]         # Note: 1
pmatrix <- scale(protein[,vars.to.use])    	# Note: 2
pcenter <- attr(pmatrix, "scaled:center")  	# Note: 3
pscale <- attr(pmatrix, "scaled:scale")

# Note 1:
#   Use all the columns except the first
#   (Country).

# Note 2:
#   The output of scale() is a matrix. For the
#   purposes of this chapter, you can think of a
#   matrix as a data frame with all numeric columns
#   (this isn’t strictly true, but it’s close enough).

# Note 3:
#   The scale() function annotates its output
#   with two attributes—scaled:center returns the mean
#   values of all the columns, and scaled:scale
#   returns the standard deviations. You’ll store
#   these away so you can “unscale” the data
#   later.

# Title: Hierarchical clustering

d <- dist(pmatrix, method="euclidean")     # Note: 1
pfit <- hclust(d, method="ward")         	# Note: 2
plot(pfit, labels=protein$Country)      	# Note: 3

# Note 1:
#   Create the distance matrix.

# Note 2:
#   Do the clustering.

# Note 3:
#   Plot the dendrogram.

# Title: Hierarchical clustering

d <- dist(pmatrix, method="euclidean")     # Note: 1
pfit <- hclust(d, method="ward")         	# Note: 2
plot(pfit, labels=protein$Country)      	# Note: 3

# Note 1:
#   Create the distance matrix.

# Note 2:
#   Do the clustering.

# Note 3:
#   Plot the dendrogram.

# Title: Extracting the clusters found by hclust()

groups <- cutree(pfit, k=5)

print_clusters <- function(labels, k) {               # Note: 1
  for(i in 1:k) {
    print(paste("cluster", i))
    print(protein[labels==i,c("Country","RedMeat","Fish","Fr.Veg")])
  }
}

print_clusters(groups, 5)

# Note 1:
#   A convenience function for printing out the
#   countries in each cluster, along with the values
#   for red meat, fish, and fruit/vegetable
#   consumption. We’ll use this function throughout
#   this section. Note that the function is hardcoded
#   for the protein dataset.

# Title: Projecting the clusters on the first two principal components

library(ggplot2)
princ <- prcomp(pmatrix)      # Note: 1
nComp <- 2
project <- predict(princ, newdata=pmatrix)[,1:nComp]      	# Note: 2
project.plus <- cbind(as.data.frame(project),             	# Note: 3
                      cluster=as.factor(groups),
                      country=protein$Country)
ggplot(project.plus, aes(x=PC1, y=PC2)) +                	# Note: 4
  geom_point(aes(shape=cluster)) +
  geom_text(aes(label=country),
            hjust=0, vjust=1)

# Note 1:
#   Calculate the principal components of the
#   data.

# Note 2:
#   The predict() function will rotate the data
#   into the space described by the principal
#   components. We only want the projection on the
#   first two axes.

# Note 3:
#   Create a data frame with the transformed
#   data, along with the cluster label and country
#   label of each point.

# Note 4:
#   Plot it.

# Title: Running clusterboot() on the protein data

library(fpc)                                    # Note: 1
kbest.p<-5                                                   	# Note: 2
cboot.hclust <- clusterboot(pmatrix,clustermethod=hclustCBI, 	# Note: 3
                            method="ward", k=kbest.p)

summary(cboot.hclust$result)                              	# Note: 4


groups<-cboot.hclust$result$partition                      	# Note: 5
print_clusters(groups, kbest.p)                           	# Note: 6


# Note 1:
#   Load the fpc package. You may have to
#   install it first. We’ll discuss installing R
#   packages in appendix .

# Note 2:
#   Set the desired number of clusters.

# Note 3:
#   Run clusterboot() with hclust
#   ('clustermethod=hclustCBI') using Ward’s method
#   ('method="ward"') and kbest.p clusters
#   ('k=kbest.p'). Return the results in an object
#   called cboot.hclust.

# Note 4:
#   The results of the clustering are in
#   cboot.hclust$result. The output of the hclust()
#   function is in cboot.hclust$result$result.

# Note 5:
#   cboot.hclust$result$partition returns a
#   vector of clusterlabels.

# Note 6:
#   The clusters are the same as those produced
#   by a direct call to hclust().

# Note 7:
#   The vector of cluster stabilities.

# Note 8:
#   The count of how many times each cluster was
#   dissolved. By default clusterboot() runs 100
#   bootstrap iterations.

# Title: Calculating total within sum of squares

sqr_edist <- function(x, y) {               # Note: 1
  sum((x-y)^2)
}

wss.cluster <- function(clustermat) {     	# Note: 2
  c0 <- apply(clustermat, 2, FUN=mean)    	# Note: 3
  sum(apply(clustermat, 1, FUN=function(row){sqr_edist(row,c0)}))     	# Note: 4
}

wss.total <- function(dmatrix, labels) {                               	# Note: 5
  wsstot <- 0
  k <- length(unique(labels))
  for(i in 1:k)
    wsstot <- wsstot + wss.cluster(subset(dmatrix, labels==i))         	# Note: 6
  wsstot
}

# Note 1:
#   Function to calculate squared distance
#   between two vectors.

# Note 2:
#   Function to calculate the WSS for a single
#   cluster, which is represented as a matrix (one row
#   for every point).

# Note 3:
#   Calculate the centroid of the cluster (the
#   mean of all the points).

# Note 4:
#   Calculate the squared difference of every
#   point in the cluster from the centroid, and sum
#   all the distances.

# Note 5:
#   Function to compute the total WSS from a set
#   of data points and cluster labels.

# Note 6:
#   Extract each cluster, calculate the
#   cluster’s WSS, and sum all the values.

# Title: The Calinski-Harabasz index

totss <- function(dmatrix) {                   # Note: 1
  grandmean <- apply(dmatrix, 2, FUN=mean)
  sum(apply(dmatrix, 1, FUN=function(row){sqr_edist(row, grandmean)}))
}


ch_criterion <- function(dmatrix, kmax, method="kmeans") {     	# Note: 2
  if(!(method %in% c("kmeans", "hclust"))) {
    stop("method must be one of c('kmeans', 'hclust')")
  }
  npts <- dim(dmatrix)[1]  # number of rows.

  totss <- totss(dmatrix)                                       	# Note: 3

  wss <- numeric(kmax)
  crit <- numeric(kmax)
  wss[1] <- (npts-1)*sum(apply(dmatrix, 2, var))                	# Note: 4
  for(k in 2:kmax) {                                           	# Note: 5
    if(method=="kmeans") {
      clustering<-kmeans(dmatrix, k, nstart=10, iter.max=100)
      wss[k] <- clustering$tot.withinss
    }else {  # hclust                                          	# Note: 6
      d <- dist(dmatrix, method="euclidean")
      pfit <- hclust(d, method="ward")
      labels <- cutree(pfit, k=k)
      wss[k] <- wss.total(dmatrix, labels)
    }
  }
  bss <- totss - wss                                            	# Note: 7
  crit.num <- bss/(0:(kmax-1))                                  	# Note: 8
  crit.denom <- wss/(npts - 1:kmax)                             	# Note: 9
  list(crit = crit.num/crit.denom, wss = wss, totss = totss)   	# Note: 10
}

# Note 1:
#   Convenience function to calculate the total
#   sum of squares.

# Note 2:
#   A function to calculate the CH index for a
#   number of clusters from 1 to kmax.

# Note 3:
#   The total sum of squares is independent of
#   the clustering.

# Note 4:
#   Calculate WSS for k=1 (which is really just
#   total sum of squares).

# Note 5:
#   Calculate WSS for k from 2 to kmax. kmeans()
#   returns the total WSS as one of its
#   outputs.

# Note 6:
#   For hclust(), calculate total WSS by
#   hand.

# Note 7:
#   Calculate BSS for k from 1 to kmax.

# Note 8:
#   Normalize BSS by k-1.

# Note 9:
#   Normalize WSS by npts - k.

# Note 10:
#   Return a vector of CH indices and of WSS for
#   k from 1 to kmax. Also return total sum of
#   squares.

# Title: Evaluating clusterings with different numbers of clusters

library(reshape2)                                           # Note: 1
clustcrit <- ch_criterion(pmatrix, 10, method="hclust")     	# Note: 2
critframe <- data.frame(k=1:10, ch=scale(clustcrit$crit),   	# Note: 3
                        wss=scale(clustcrit$wss))
critframe <- melt(critframe, id.vars=c("k"),                	# Note: 4
                  variable.name="measure",
                  value.name="score")
ggplot(critframe, aes(x=k, y=score, color=measure)) +     	# Note: 5
  geom_point(aes(shape=measure)) + geom_line(aes(linetype=measure)) +
  scale_x_continuous(breaks=1:10, labels=1:10)

# Note 1:
#   Load the reshape2 package (for the melt()
#   function).

# Note 2:
#   Calculate both criteria for 1–10
#   clusters.

# Note 3:
#   Create a data frame with the number of
#   clusters, the CH criterion, and the WSS criterion.
#   We’ll scale both the CH and WSS criteria to
#   similar ranges so that we can plot them both on
#   the same graph.

# Note 4:
#   Use the melt() function to put the data
#   frame in a shape suitable for ggplot

# Note 5:
#   Plot it.

# Title: Running k-means with k=5

pclusters <- kmeans(pmatrix, kbest.p, nstart=100, iter.max=100)     # Note: 1
summary(pclusters)                                               	# Note: 2
pclusters$centers                                                	# Note: 3
pclusters$size                                                  	# Note: 4
groups <- pclusters$cluster                                      	# Note: 5
print_clusters(groups, kbest.p)                                 	# Note: 6

# Note 1:
#   Run kmeans() with five clusters (kbest.p=5),
#   100 random starts, and 100 maximum iterations per
#   run.

# Note 2:
#   kmeans() returns all the sum of squares
#   measures.

# Note 3:
#   pclusters$centers is a matrix whose rows are
#   the centroids of the clusters. Note that
#   pclusters$centers is in the scaled coordinates,
#   not the original protein coordinates.

# Note 4:
#   pclusters$size returns the number of points
#   in each cluster. Generally (though not always) a
#   good clustering will be fairly well balanced: no
#   extremely small clusters and no extremely large
#   ones.

# Note 5:
#   pclusters$cluster is a vector of cluster
#   labels.

# Note 6:
#   In this case, kmeans() and hclust() returned
#   the same clustering. This won’t always be
#   true.

# Title: Plotting cluster criteria

clustering.ch <- kmeansruns(pmatrix, krange=1:10, criterion="ch")     # Note: 1
clustering.ch$bestk                                                	# Note: 2

clustering.asw <- kmeansruns(pmatrix, krange=1:10, criterion="asw") 	# Note: 3
clustering.asw$bestk

clustering.ch$crit                                                 	# Note: 4

clustcrit$crit                                                     	# Note: 5

critframe <- data.frame(k=1:10, ch=scale(clustering.ch$crit),        	# Note: 6
                          asw=scale(clustering.asw$crit))
critframe <- melt(critframe, id.vars=c("k"),
                    variable.name="measure",
                    value.name="score")
ggplot(critframe, aes(x=k, y=score, color=measure)) +
  geom_point(aes(shape=measure)) + geom_line(aes(linetype=measure)) +
  scale_x_continuous(breaks=1:10, labels=1:10)
summary(clustering.ch)                                            	# Note: 7

# Note 1:
#   Run kmeansruns() from 1–10 clusters, and the
#   CH criterion. By default, kmeansruns() uses 100
#   random starts and 100 maximum iterations per
#   run.

# Note 2:
#   The CH criterion picks two clusters.

# Note 3:
#   Run kmeansruns() from 1–10 clusters, and the
#   average silhouette width criterion. Average
#   silhouette width picks 3 clusters.

# Note 4:
#   The vector of criterion values is called
#   crit.

# Note 5:
#   Compare the CH values for kmeans() and
#   hclust(). They’re not quite the same, because the
#   two algorithms didn’t pick the same
#   clusters.

# Note 6:
#   Plot the values for the two criteria.

# Note 7:
#   kmeansruns() also returns the output of
#   kmeans for k=bestk.

# Title: Running clusterboot() with k-means

kbest.p<-5
cboot<-clusterboot(pmatrix, clustermethod=kmeansCBI,
                   runs=100,iter.max=100,
                   krange=kbest.p, seed=15555)                 # Note: 1

groups <- cboot$result$partition
print_clusters(cboot$result$partition, kbest.p)

cboot$bootmean

cboot$bootbrd


# Note 1:
#   We’ve set the seed for the random generator
#   so the results are reproducible.

# Title: A function to assign points to a cluster

assign_cluster <- function(newpt, centers, xcenter=0, xscale=1) {   # Note: 1
  xpt <- (newpt - xcenter)/xscale                                	# Note: 2
  dists <- apply(centers, 1, FUN=function(c0){sqr_edist(c0, xpt)})  	# Note: 3
  which.min(dists)                                                 	# Note: 4
}

# Note 1:
#   A function to assign a new data point newpt to
#   a clustering described by centers, a matrix where
#   each row is a cluster centroid. If the data was
#   scaled (using scale()) before clustering, then
#   xcenter and xscale are the scaled:center and
#   scaled:scale attributes, respectively.

# Note 2:
#   Center and scale the new data point.

# Note 3:
#   Calculate how far the new data point is from
#   each of the cluster centers.

# Note 4:
#   Return the cluster number of the closest
#   centroid.

# Title: An example of assigning points to cluster

rnorm.multidim <- function(n, mean, sd, colstr="x") {      # Note: 1
  ndim <- length(mean)
  data <- NULL
  for(i in 1:ndim) {
    col <- rnorm(n, mean=mean[[i]], sd=sd[[i]])
    data<-cbind(data, col)
  }
  cnames <- paste(colstr, 1:ndim, sep='')
  colnames(data) <- cnames
  data
}

mean1 <- c(1, 1, 1)                    	# Note: 2
sd1 <- c(1, 2, 1)

mean2 <- c(10, -3, 5)
sd2 <- c(2, 1, 2)

mean3 <- c(-5, -5, -5)
sd3 <- c(1.5, 2, 1)

clust1 <- rnorm.multidim(100, mean1, sd1)           	# Note: 3
clust2 <- rnorm.multidim(100, mean2, sd2)
clust3 <- rnorm.multidim(100, mean3, sd3)
toydata <- rbind(clust3, rbind(clust1, clust2))

tmatrix <- scale(toydata)                          	# Note: 4
tcenter <- attr(tmatrix, "scaled:center")        	# Note: 5
tscale<-attr(tmatrix, "scaled:scale")
kbest.t <- 3
tclusters <- kmeans(tmatrix, kbest.t, nstart=100, iter.max=100)   	# Note: 6

tclusters$size              	# Note: 7


unscale <- function(scaledpt, centervec, scalevec) {    	# Note: 8
  scaledpt*scalevec + centervec
}

unscale(tclusters$centers[1,], tcenter, tscale)   	# Note: 9

mean2

unscale(tclusters$centers[2,], tcenter, tscale)   	# Note: 10

mean3

unscale(tclusters$centers[3,], tcenter, tscale)    	# Note: 11

mean1


assign_cluster(rnorm.multidim(1, mean1, sd1),  	# Note: 12
                 tclusters$centers,
                 tcenter, tscale)
                                                	# Note: 13


assign_cluster(rnorm.multidim(1, mean2, sd1),   	# Note: 14
                 tclusters$centers,
                 tcenter, tscale)
                                              	# Note: 15


assign_cluster(rnorm.multidim(1, mean3, sd1),     	# Note: 16
                 tclusters$centers,
                 tcenter, tscale)
                                          	# Note: 17


# Note 1:
#   A function to generate n points drawn from a
#   multidimensional Gaussian distribution with
#   centroid mean and standard deviation sd. The
#   dimension of the distribution is given by the
#   length of the vector mean.

# Note 2:
#   The parameters for three Gaussian
#   distributions.

# Note 3:
#   Create a dataset with 100 points each drawn
#   from the above distributions.

# Note 4:
#   Scale the dataset.

# Note 5:
#   Store the centering and scaling parameters for
#   future use.

# Note 6:
#   Cluster the dataset, using k-means with three
#   clusters.

# Note 7:
#   The resulting clusters are about the right
#   size.

# Note 8:
#   A function to “unscale” data points (put them
#   back in the coordinates of the original
#   dataset).

# Note 9:
#   Unscale the first centroid. It corresponds to
#   our original distribution 2.

# Note 10:
#   The second centroid corresponds to the
#   original distribution 3.

# Note 11:
#   The third centroid corresponds to the original
#   distribution 1.

# Note 12:
#   Generate a random point from the original
#   distribution 1 and assign it to one of the
#   discovered clusters.

# Note 13:
#   It’s assigned to cluster 3, as we would
#   expect.

# Note 14:
#   Generate a random point from the original
#   distribution 2 and assign it.

# Note 15:
#   It’s assigned to cluster 1.

# Note 16:
#   Generate a random point from the original
#   distribution 3 and assign it.

# Note 17:
#   It’s assigned to cluster 2.

##############################
# ASSOCIATION RULES EXAMPLE: #
##############################

# Title: Reading in the book data
setwd("~/Desktop/Data_Science/zmPDSwR-master/Bookdata")
library(arules)    # Note: 1
bookbaskets <- read.transactions("bookdata.tsv.gz", format="single",  	# Note: 2
                                 sep="\t",                    	# Note: 3
                                 cols=c("userid", "title"),    	# Note: 4
                                 rm.duplicates=T)       	# Note: 5

# Note 1:
#   Load the arules package.

# Note 2:
#   Specify the file and the file format.

# Note 3:
#   Specify the column separator (a tab).

# Note 4:
#   Specify the column of transaction IDs and of
#   item IDs, respectively.

# Note 5:
#   Tell the function to look for and remove
#   duplicate entries (for example, multiple entries
#   for “The Hobbit” by the same user).

# Title: Examining the transaction data

class(bookbaskets)               # Note: 1

bookbaskets                    	# Note: 2

dim(bookbaskets)               	# Note: 3

colnames(bookbaskets)[1:5]     	# Note: 4

rownames(bookbaskets)[1:5]        	# Note: 5


# Note 1:
#   The object is of class transactions.

# Note 2:
#   Printing the object tells you its
#   dimensions.

# Note 3:
#   You can also use dim() to see the dimensions
#   of the matrix.

# Note 4:
#   The columns are labeled by book
#   title.

# Note 5:
#   The rows are labeled by customer.


# Title: Examining the size distribution

quantile(basketSizes, probs=seq(0,1,0.1))       # Note: 1

 library(ggplot2)                              	# Note: 2
 ggplot(data.frame(count=basketSizes)) +
  geom_density(aes(x=count), binwidth=1) +
  scale_x_log10()

# Note 1:
#   Look at the basket size distribution, in 10%
#   increments.

# Note 2:
#   Plot the distribution to get a better
#   look.

# (informalexample 8.8 of section 8.2.3)  : Unsupervised methods : Association rules : Mining association rules with the arules package

 bookFreq <- itemFrequency(bookbaskets)
summary(bookFreq)

 sum(bookFreq)

# Title: Finding the ten most frequent books

 bookCount <- (bookFreq/sum(bookFreq))*sum(basketSizes)       # Note: 1
 summary(bookCount)

 orderedBooks <- sort(bookCount, decreasing=T)   	# Note: 2
 orderedBooks[1:10]

 orderedBooks[1]/dim(bookbaskets)[1]                 	# Note: 3


# Note 1:
#   Get the absolute count of book
#   occurrences.

# Note 2:
#   Sort the count and list the 10 most popular
#   books.

# Note 3:
#   The most popular book in the dataset
#   occurred in fewer than 3% of the baskets.


# Title: Finding the association rules

 rules <- apriori(bookbaskets_use,                                    # Note: 1
                   parameter =list(support = 0.002, confidence=0.75))

 summary(rules)
# set of 191 rules                            	# Note: 2
#
# rule length distribution (lhs + rhs):sizes         	# Note: 3
# 2   3   4   5
# 11 100  66  14
#
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 2.000   3.000   3.000   3.435   4.000   5.000

# summary of quality measures:                           	# Note: 4
#   support           confidence          lift
# Min.   :0.002009   Min.   :0.7500   Min.   : 40.89
# 1st Qu.:0.002131   1st Qu.:0.8113   1st Qu.: 86.44
# Median :0.002278   Median :0.8468   Median :131.36
# Mean   :0.002593   Mean   :0.8569   Mean   :129.68
# 3rd Qu.:0.002695   3rd Qu.:0.9065   3rd Qu.:158.77
# Max.   :0.005830   Max.   :0.9882   Max.   :321.89

#mining info:                                           	# Note: 5
#  data ntransactions support confidence
#bookbaskets_use         40822   0.002       0.75

# Note 1:
#   Call apriori() with a minimum support of
#   0.002 and a minimum confidence of 0.75.

# Note 2:
#   The summary of the apriori() output reports
#   the number of rules found;...

# Note 3:
#   ...the distribution of rule lengths (in this
#   example, most rules contain 3 items—2 on the left
#   side, X (lhs), and one on the right side, Y
#   (rhs));...

# Note 4:
#   ...a summary of rule quality measures,
#   including support and confidence;...

# Note 5:
#   ...and some information on how apriori() was
#   called.

# Title: Scoring rules

measures <- interestMeasure(rules,                              # Note: 1
        method=c("coverage", "fishersExactTest"),    	# Note: 2
        transactions=bookbaskets_use)                	# Note: 3
summary(measures)
# coverage        fishersExactTest
# Min.   :0.002082   Min.   : 0.000e+00
# 1st Qu.:0.002511   1st Qu.: 0.000e+00
# Median :0.002719   Median : 0.000e+00
# Mean   :0.003039   Mean   :5.080e-138
# 3rd Qu.:0.003160   3rd Qu.: 0.000e+00
# Max.   :0.006982   Max.   :9.702e-136

# Note 1:
#   The call to interestMeasure() takes as
#   arguments the discovered rules,...

# Note 2:
#   ...a list of interest measures to
#   apply,...

# Note 3:
#   ...and a dataset to evaluate the interest
#   measures over. This is usually the same set used
#   to mine the rules, but it needn’t be. For
#   instance, you can evaluate the rules over the full
#   dataset, bookbaskets, to get coverage estimates
#   that reflect all the customers, not just the ones
#   who showed interest in more than one book.

# Title: Finding rules with restrictions

brules <- apriori(bookbaskets_use,
                  parameter =list(support = 0.001,      # Note: 1
                                  confidence=0.6),
                  appearance=list(rhs=c("The Lovely Bones: A Novel"),  	# Note: 2
                                  default="lhs"))                      	# Note: 3
 summary(brules)
# set of 46 rules
#
# rule length distribution (lhs + rhs):sizes
# 3  4
# 44  2
#
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 3.000   3.000   3.000   3.043   3.000   4.000
#
# summary of quality measures:
#   support           confidence          lift
# Min.   :0.001004   Min.   :0.6000   Min.   :21.81
# 1st Qu.:0.001029   1st Qu.:0.6118   1st Qu.:22.24
# Median :0.001102   Median :0.6258   Median :22.75
# Mean   :0.001132   Mean   :0.6365   Mean   :23.14
# 3rd Qu.:0.001219   3rd Qu.:0.6457   3rd Qu.:23.47
# Max.   :0.001396   Max.   :0.7455   Max.   :27.10
#
# mining info:
#   data ntransactions support confidence
# bookbaskets_use         40822   0.001        0.6

# Note 1:
#   Relax the minimum support to 0.001 and the
#   minimum confidence to 0.6.

# Note 2:
#   Only The Lovely Bones
#   is allowed to appear on the right side of the
#   rules.

# Note 3:
#   By default, all the books can go into the
#   left side of the rules.

# Title: Inspecting rules

brulesConf <- sort(brules, by="confidence")    # Note: 1

inspect(head(lhs(brulesConf), n=5))      	# Note: 2
# items
# 1 {Divine Secrets of the Ya-Ya Sisterhood: A Novel,
#    Lucky : A Memoir}
# 2 {Lucky : A Memoir,
#    The Notebook}
# 3 {Lucky : A Memoir,
#    Wild Animus}
# 4 {Midwives: A Novel,
#    Wicked: The Life and Times of the Wicked Witch of the West}
# 5 {Lucky : A Memoir,
#    Summer Sisters}

# Note 1:
#   Sort the rules by confidence.

# Note 2:
#   Use the lhs() function to get the left
#   itemsets of each rule; then inspect the top
#   five.

# Title: Inspecting rules with restrictions

brulesSub <- subset(brules, subset=!(lhs %in% "Lucky : A Memoir"))    # Note: 1
brulesConf <- sort(brulesSub, by="confidence")

 inspect(head(lhs(brulesConf), n=5))
# items
# 1 {Midwives: A Novel,
#    Wicked: The Life and Times of the Wicked Witch of the West}
# 2 {She's Come Undone,
#    The Secret Life of Bees,
#    Wild Animus}
# 3 {A Walk to Remember,
#    The Nanny Diaries: A Novel}
# 4 {Beloved,
#    The Red Tent}
# 5 {The Da Vinci Code,
#    The Reader}

# Note 1:
#   Restrict to the subset of rules where
#   Lucky is not in the left
#   side.
