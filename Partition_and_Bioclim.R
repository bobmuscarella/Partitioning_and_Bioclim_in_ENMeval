######################################
####### DATA PARTITIONING STEP #######
######################################
# The 'get.group.data' function below depends on ENMeval functions to partition data for model evaluation.
# Currently in Wallace, background localities are generated when ENMeval is run but this needs to be moved forward to enable other algorithms to work.
# See documentation of ENMeval functions for more information.

library(ENMeval)

get.group.data <- function(method, occ, env, bg.coords=NULL, aggregation.factor=NULL, kfolds=NULL, occ.grp=NULL, bg.grp=NULL){
	if(method != 'user') { occ <- as.data.frame(occ) }
	if(is.null(bg.coords)) { 
		bg.coords <- randomPoints(env, 10000)
		}
	bg.coords <- as.data.frame(bg.coords)
	if(method=='block') { group.data <- get.block(occ, bg.coords) }
	if(method=='checkerboard1') { group.data <- get.checkerboard1(occ, env, bg.coords, aggregation.factor) }
	if(method=='checkerboard2') { group.data <- get.checkerboard2(occ, env, bg.coords, aggregation.factor) }
	if(method=='jackknife') { group.data <- get.jackknife(occ, bg.coords) }
	if(method=='randomkfold') { group.data <- get.randomkfold(occ, bg.coords, kfolds) }
	if(method=='user') { group.data <- get.user(occ.grp, bg.grp) }
	return(list(occ.pts=occ, bg.pts=bg.coords, occ.grp=group.data[[1]], bg.grp=group.data[[2]]))
}


#############################################
####### TEST THE PARTITIONING METHODS ######
#############################################

### SIMULATE SOME DATA ###
r1 <- raster(matrix(nrow=50, ncol=50, data=runif(10000, 0, 25)))
r2 <- raster(matrix(nrow=50, ncol=50, data=rep(1:100, each=100), byrow=TRUE))
r3 <- raster(matrix(nrow=50, ncol=50, data=rep(1:100, each=100)))
r4 <- raster(matrix(nrow=50, ncol=50, data=c(rep(1,1000),rep(2,500)),byrow=TRUE))
values(r4) <- as.factor(values(r4))
env <- stack(r1,r2,r3,r4)
nocc <- 50
x <- (rpois(nocc, 2) + abs(rnorm(nocc)))/11
y <- runif(nocc, 0, .99)
occ <- cbind(x,y)


### RUN PARTITIONING METHODS ###
# BLOCK
group.data <- get.group.data(method='block', env=env, occ=occ, bg.coords=bg.coords)

# CHECKERBOARD1
group.data <- get.group.data(method='checkerboard1', occ=occ, env=env, bg.coords=bg.coords, aggregation.factor=2)

# CHECKERBOARD2
group.data <- get.group.data(method='checkerboard2', occ=occ, env=env, bg.coords=bg.coords, aggregation.factor=2)

# JACKKNIFE
group.data <- get.group.data(method='jackknife', occ=occ, env=env, bg.coords=bg.coords)

# RANDOM K FOLD
group.data <- get.group.data(method='randomkfold', occ=occ, env=env, bg.coords=bg.coords, kfolds=5)

# USER
occ.grp=sample(1:4, nrow(occ), replace=T)
bg.grp=sample(1:4, nrow(bg.coords), replace=T)
group.data <- get.group.data(method='user', occ=occ, env=env, occ.grp=occ.grp, bg.grp=bg.grp)



########################################
####### RUN AND EVALUATE BIOCLIM ######
########################################
# This uses the dismo implementation of bioclim models and the partitioned data from above.
# The bioclim model does not require background points but they are used to generate some of the evaluation statistics (AUC, AUC.diff)


BioClim_eval <- function (group.data, env) {

# RUN FULL DATA MODEL
	full.mod <- bioclim(env, group.data$occ.pts)
	pred <- predict(env, full.mod)

# CREATE HOLDERS FOR RESULTS
	AUC.TEST <- double()
	AUC.DIFF <- double()
	OR10 <- double()
	ORmin <- double()

# SET NUMBER OF TEST BINS
	nk <- length(unique(group.data$occ.grp))

	for (k in 1:nk) {

	# SPLIT TEST AND TRAIN DATA
		test.pts <- occ[group.data$occ.grp == k, ]
		train.pts <- occ[group.data$occ.grp != k, ]
		bg.pts <- group.data$bg.pts[group.data$bg.grp != k, ]
		mod <- bioclim(env, train.pts)
		
	# GET AUC METRICS
		AUC.TEST[k] <- evaluate(p=test.pts, a=bg.pts, mod=mod, x=env)@auc
		AUC.TRAIN <- evaluate(p=train.pts, a=bg.pts, mod=mod, x=env)@auc
		AUC.DIFF[k] <- max(0, AUC.TRAIN - AUC.TEST[k])
	
	# GET PREDICTED VALUES AT OCCURRENCES FOR OMISSION RATE STATS
		train.pred <- predict(env, mod)
		p.train <- extract(train.pred, train.pts)
		p.test <- extract(train.pred, test.pts)
	
	# FIND THRESHOLD FOR OR10
		if (nrow(train.pts) < 10) {
		n90 <- floor(nrow(train.pts) * 0.9)
		} else {
		n90 <- ceiling(nrow(train.pts) * 0.9)
		}

	# GET OMISSION RATE STATS
		train.thr.10 <- rev(sort(p.train))[n90]
		OR10[k] <- mean(p.test < train.thr.10)
		ORmin[k] <- mean(p.test < min(p.train))
	}

# COMPILE AND SUMMARIZE RESULTS
stats <- as.data.frame(rbind(AUC.DIFF, AUC.TEST, OR10, ORmin))
stats <- cbind(apply(stats, 1, mean), corrected.var(stats, nk), stats)
colnames(stats) <- c("Mean", "Variance", paste("Bin", 1:nk))
rownames(stats) <- c("AUC.DIFF", "AUC.TEST","OR10","ORmin")

# THIS FORMAT FOR RETURNED DATA ATTEMPTS TO MATCH WHAT HAPPENS IN WALLACE ALREADY FOR ENMEVAL.
return(list(full.mod, stats, pred))  
}

###########################
####### TRY IT OUT! #######
###########################

results <- BioClim_eval(group.data, env)

results

plot(results[[3]]); points(group.data$occ.pts)
