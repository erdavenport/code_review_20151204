# README

## Contents:
This repo contains the script that we covered in code review on 12-4-2015: `beta_diversity_permutations_one_twin_120115ERD.R`. 
In addition:  
* `scripts/`:
  * `beta_diversity_permutations_one_twin_120115ERD.R`: The original script (with some notes in it from the CR)  
  * `beta_diversity_permutations_one_twin_UPDATED_120115ERD.R`: The new script with comments incorporated   
  * `submit_beta_diversity_permutations_one_twin_1200115.sh`: The original submission script to run the permutations on the cluster  
  * `submit_beta_diversity_permutations_one_twin_UPDATED_1200115.sh`: The new submission script to run the UPDATED script on the cluster  
  * `code_review_notes_to_Emily.r`: File Elissa sent with her updates to the script. Elissa == awesome.  
* `results/`:  
  * This folder contains all of the data files necessary for running the script and will hold the output of the permutations.  
* `code_review_setup.pptx`:  
  * The powerpoint I ran through quickly at the beginning of code review to set up what the script was doing.  

## To run:
To run the original submission script, run the following command from the `scripts/` directory:  

```
for i in pABO
do
qsub submit_beta_diversity_permutations_one_twin_1200115.sh \
	--abo_file=../results/$i/abo.table.$i.txt \
	--cov_file=../results/$i/covs.table.$i.txt \
	--bc_file=../results/$i/beta.div.bc.$i.txt \
	--uu_file=../results/$i/beta.div.wu.$i.txt \
	--wu_file=../results/$i/beta.div.uu.$i.txt \
	--n=1000 \
	--nproc=8 \
	--outpath=../results/$i/8_beta_diversity/
done
```

To run the new submission script, run the following command from the `scripts/` directory:  

```
for i in pABO
do
qsub submit_beta_diversity_permutations_one_twin_UPDATED_1200115.sh \
	--abo_file=../results/$i/abo.table.$i.txt \
	--cov_file=../results/$i/covs.table.$i.txt \
	--bc_file=../results/$i/beta.div.bc.$i.txt \
	--uu_file=../results/$i/beta.div.wu.$i.txt \
	--wu_file=../results/$i/beta.div.uu.$i.txt \
	--n=1000 \
	--nproc=8 \
	--outpath=../results/$i/8_beta_diversity/
done
```

## Notes:  
* Create one melted look-up table that contains all of the distance metrics, rather than having three tables floating around:  
Original (and UGLY):  
```
bc2[upper.tri(bc2, diag=TRUE)] <- 999
colnames(bc2) <- abo2$ABO
rownames(bc2) <- abo2$ABO
bc.abo.M <- subset(melt(bc2), value!=999)
colnames(bc2) <- abo2$secretor_status
rownames(bc2) <- abo2$secretor_status
bc.ss.M <- subset(melt(bc2), value!=999)

wu2[upper.tri(wu2, diag=TRUE)] <- 999
colnames(wu2) <- abo2$ABO
rownames(wu2) <- abo2$ABO
wu.abo.M <- subset(melt(wu2), value!=999)
colnames(wu2) <- abo2$secretor_status
rownames(wu2) <- abo2$secretor_status
wu.ss.M <- subset(melt(wu2), value!=999)

uu2[upper.tri(uu2, diag=TRUE)] <- 999
colnames(uu2) <- abo2$ABO
rownames(uu2) <- abo2$ABO
uu.abo.M <- subset(melt(uu2), value!=999)
colnames(uu2) <- abo2$secretor_status
rownames(uu2) <- abo2$secretor_status
uu.ss.M <- subset(melt(uu2), value!=999)
```

Updated:  
```
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
```

* Replace the function to determine "within" and "between" with ifelse function (hat tip to Elissa, I didn't know about this!):  

Original:  
```
groupme <- function(m) {
	if (is.na(m[1]) | is.na(m[2])) { # for missing SS individuals
		return(NA)
	} else if (m[1] == m[2]) {
		return("within")
	} else {
		return("between")
	}
}

abo.group <- apply(bc.abo.M, 1, groupme)
ss.group <- apply(bc.ss.M, 1, groupme)
```

Updated:  
```
abo.group <- ifelse(abo.i1 == abo.i2, "within", "btwn")
ss.group <- ifelse(ss.i1 == ss.i2, "within", "btwn")
```

* Instead of melting the distance matrix into a three column table every permutation, melt once and then do a look up to match ID names:

Updated:  
```
# generate sample permutations
ID.idx <- 1:ncol(bc2)
set.seed(1)
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
```

Big big thanks to everyone for comments!