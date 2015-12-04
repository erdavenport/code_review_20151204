
## EJC: just adding one library here

##### Load libraries:
library("testit")
suppressMessages(library("foreach"))
suppressMessages(library("doParallel"))
library("reshape2")
#ADD:
library(parallel)

#~~~~
#~~~~

## EJC: corresponding to lines 110-217

##### Create melted distance table for each distance metric with either ABO or SS labels:

bc2[upper.tri(bc2, diag=T)] <- NA
wu2[upper.tri(wu2, diag=T)] <- NA
uu2[upper.tri(uu2, diag=T)] <- NA

melt.dist <- melt(bc2)
names(melt.dist) <- c("ID1", "ID2", "bc2")

# add other distance metrics
melt.dist$wu2 <- melt(wu2)$value
melt.dist$uu2 <- melt(uu2)$value

# remove NA rows
melt.dist <- melt.dist[!is.na(melt.dist$bc2), ]

# map ABO and secretor status
i1.idx <- match(melt.dist$ID1, colnames(bc2))
i2.idx <- match(melt.dist$ID2, colnames(bc2))
abo.i1 <- abo2$ABO[i1.idx]
abo.i2 <- abo2$ABO[i2.idx]
ss.i1 <- abo2$secretor_status[i1.idx]
ss.i2 <- abo2$secretor_status[i2.idx]

# make group calls
abo.group <- ifelse(abo.i1 == abo.i2, "within", "btwn")
ss.group <- ifelse(ss.i1 == ss.i2, "within", "btwn")


##### Calculate actual t-test p-value:

abo.p <- sapply(melt.dist[, 3:5], function(d) t.test(d ~ abo.group)$p.value)
ss.p <- sapply(melt.dist[, 3:5], function(d) t.test(d ~ ss.group)$p.value)

# NOTE: different order of p-values
actual.ps <- c(abo.p, ss.p)
names(actual.ps) <- c("bc_abo", "wu_abo", "uu_abo", "bc_ss", "wu_ss", "uu_ss")

##### Run the permutations:
print("running permutations in parallel")

# Start time:
strt <- Sys.time()

# generate sample permutations
ID.idx <- 1:ncol(bc2)
perm.list <- lapply(1:n, function(x) sample(ID.idx))

perm.out <- mclapply(perm.list, function(j)
{	
	# Apply permuted sample indices
	i1.idx <- match(melt.dist$ID1, colnames(bc2)[j])
	i2.idx <- match(melt.dist$ID2, colnames(bc2)[j])

	# map ABO and secretor status
	abo.i1 <- abo2$ABO[i1.idx]
	abo.i2 <- abo2$ABO[i2.idx]
	ss.i1 <- abo2$secretor_status[i1.idx]
	ss.i2 <- abo2$secretor_status[i2.idx]

	# make group calls
	abo.group <- ifelse(abo.i1 == abo.i2, "within", "btwn")
	ss.group <- ifelse(ss.i1 == ss.i2, "within", "btwn")

	# Calculate permutation p-values for each distance metric and ABO/SS:
	abo.p <- sapply(melt.dist[, 3:5], function(d) t.test(d ~ abo.group)$p.value)
	ss.p <- sapply(melt.dist[, 3:5], function(d) t.test(d ~ ss.group)$p.value)
	
	return(c(abo.p, ss.p))
	
}, mc.cores=nproc)

print(Sys.time() - strt)
stopCluster(cl)


##### Calculate permutation p-value:

print("calculating permutation p-values")

# NOTE: different order of p-values
permutations <- do.call(rbind, perm.out)
colnames(permutations) <- c("bc_abo", "wu_abo", "uu_abo", "bc_ss", "wu_ss", "uu_ss")

permutation.p <- sapply(1:length(actual.ps), function(ii) length(which(permutations[, ii] <= actual.ps[ii]))/n)

all.ps <- rbind(permutation.p, actual.ps, permutations)
rownames(all.ps)[1:2] <- c("permutation_p", "t.test_p")

#~~~~
#~~~~

