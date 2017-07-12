### ------------------------------------------------------------------------------------------------------------ ###
### Author: Anna F. Bush                                                                                         ###
### Date last revised: 15th July, 2016                                                                           ###
### ------------------------------------------------------------------------------------------------------------ ###
### 
###                                 OBSERVED MODEL SHOWING REALISTIC DISEASE DYNAMICS
###
### ------------------------------------------------------------------------------------------------------------ ###
### A program introducing a more realistic representation of disease dynamics, than expected model.
### A binary infected/uninfected system is used, where a single and randomly chosen infected individual 
### hereon named the 'leper' individual, is selected from the population, and disease spread is followed over time.
### ------------------------------------------------------------------------------------------------------------ ###
###                                              IMPORT PACKAGES                                                 ###
### ------------------------------------------------------------------------------------------------------------ ###

library(spatstat)           # spatial point patterns
library(stats)              # rbinom, rpois
library(stpp)               # rpcp
library(igraph)             # functions to plot network
library(sna)                # network metrics
library(Matrix)             # create rectangular matrix object
library(lme4)               # mixed-effects model
library(qgraph)             # network metrics          
library(tnet)               # visualise weighted network on separate screen
library(lubridate)          # time
library(binhf)              # binomial functions
library(CePa)               # network metrics

### ------------------------------------------------------------------------------------------------------------ ###
###                                                SET SYSTEM                                                    ###
### ------------------------------------------------------------------------------------------------------------ ###

setwd("E:/Social Network Project")              ### path to save results
dispersion <- c(0.1, 0.125, 0.15, 0.2, 0.3, 1)  ### list of dispersion values to run through
cluster <- "normal"                             ### normal or uniform                                 
n <- 50                                         ### number of individuals                        
p <- 1*10^-3                                    ### infectiveness of disease                     
infection.constant <- 0.5                       ### percentage of infected individuals at cut off   
n.sim <- 1                                      ### number of model iterations         
n.parents <- 10                                 ### number of parents when setting spatial field
show.spatial.field <- TRUE                      ### show spatial field and infection graph for each simulation?                                                      
show.network <- FALSE                           ### show network diagram for each simulation? (warning, slow)
show.metrics <- TRUE                            ### plot graphs of metrics against rate? 
save.results <- TRUE                            ### save metric results as csv?

### ------------------------------------------------------------------------------------------------------------ ###
###                                           DECLARE FUNCTIONS                                                  ###
### ------------------------------------------------------------------------------------------------------------ ###

### outputs an empty 3d array ready to store metrics
set.metrics <- function(dispersion, metric.names, n.sim, y.names){
  metrics <- array(NA,                                                               # NA to return error if value not filled 
                   dim=c(n.sim, length(metric.names)+3, length(dispersion)),         # n.sim x metric.names x dispersion list
                   dimnames=list(1:n.sim, c(metric.names, y.names),dispersion))      # Name the dimensions
  return(metrics)
}

### generates a random spatial field using rpcp function
### given the cluster type, dispersion value and number of individuals
### cartesian coordinates are returned
get.coords <- function(cluster, d, n, n.parents){
  coords <- rpcp(nparents=n.parents, npoints=n, lambda=0.5, nsim=1, cluster=cluster, dispersion=d, infection=FALSE)
  coords <- cbind(as.numeric(coords$xyt[,1]), as.numeric(coords$xyt[,2]))
  colnames(coords)[1:2] <- c("x", "y")                                               # rename column names to x and y
  return(coords)
}

### coords dataframe is inputted
### toroidal distances are calculated between each individual
### a euclidean distance matrix is returned
get.toroidal.distance <- function(coords){
  toroidal.distance <- matrix(0, nrow=n, ncol=n) 
  x <- coords[,1]
  y <- coords[,2]
  n <- length(x)
  for(i in 1:n){
    for(j in i:n){
      delta.x <- abs(x[i]-x[j])                        # delta.x = x distance between two coordinates (abs = absolute)
      delta.y <- abs(y[i]-y[j])                        # delta.y = y distance between two coordinates
      delta.x <- min(delta.x,1-delta.x)        
      delta.y <- min(delta.y,1-delta.y)
      dist <- sqrt(delta.x^2+delta.y^2)                # Pythagoras
      toroidal.distance[i, j] <- dist
      toroidal.distance[j, i] <- dist
    }
  }
  return(toroidal.distance)
}

### the distance matrix and number of individuals is inputted
### the relationship between distance and mean number of interactions is exponential
### the mean is inputted into the rpois function
### the expected number of interactions between i and j is returned
get.interaction.expected <- function(toroidal.distance, n){
  interaction.expected <- matrix(0, nrow=n, ncol=n)
  for (i in 1:n){
    for (j in i:n){
      if (i != j){
        interaction <- rpois(1, 4*exp(-4*toroidal.distance[i, j]))
        interaction.expected[i, j] <- interaction
        interaction.expected[j, i] <- interaction
      }
    }
  }
  return(interaction.expected)
}

### keeps track of the ongoing score of interactions between each i and j over multiple permutations
get.interaction.total <- function(interaction.expected, interaction.total){
  interaction.total <- interaction.total + interaction.expected
  return(interaction.total)
}

### inputs interaction expected, number of individuals and infectivity of disease
### multiplies the interaction.expected matrix by the infection.status matrix
### infection.status is a random deviate generate by the rbinom function
get.infection.status <- function(interaction.expected, n, p){
  total.infectious.ints<-interaction.expected%*%infection.status
  infection.status<-apply(cbind(rbinom(n,1,(1-(1-p)^total.infectious.ints)), infection.status), 1, max)
}

### calculates rate
get.rate <- function(t.max){
#  y1 <- infection.mean[t.max-1]
#  y2 <- infection.mean[t.max]
#  time <- (t.max-1)+(infection.constant-y1)/(y2-y1) # Linear interpolation
  Rate <- 1/t.max
  return(Rate)
}

### calculates lag
get.lag <- function(infection.mean){
#  y1 <- infection.mean[t.max-1]
#  y2 <- infection.mean[t.max]
#  time <- (t.max-1)+(infection.constant-y1)/(y2-y1) # Linear interpolation
  Lag <- length(infection.mean[infection.mean==1/n])
  return(Lag)
}

### calculates postlag
get.postlag <- function(infection.mean){
#  y1 <- infection.mean[t.max-1]
#  y2 <- infection.mean[t.max]
#  time <- (t.max-1)+(infection.constant-y1)/(y2-y1) # Linear interpolation
  Postlag <- 1/length(infection.mean[infection.mean>1/n])
  return(Postlag)
}

### inputs the interaction.mean and number of individuals 
### converts all xy coordinate positions into two vectors, one for x, one for y 
### the edge list is weighted by the interaction mean 
### used when plotting the network
### returns dataframe of nodes and their corresponding edge weights
get.edge.list <- function(interaction.mean, n){
  v1 <- vector(mode="numeric", length=0)
  v2 <- vector(mode="numeric", length=0)
  weight <- vector(mode="numeric", length=0)
  for(i in 1:n){
    for(j in i:n){
      if(i != j)
      {
        v1 <- append(v1, i)
        v2 <- append(v2, j)
        weight <- append(weight, interaction.mean[i, j])
      }
    }
  }
  edge.list <- as.matrix(cbind(v1, v2, weight), ncol=3)
  return(edge.list)
}

### the interaction mean is the total number of interactions 
### divided by the time in which these interactions took place
get.interaction.mean <- function(interaction.total, t.max){
  interaction.mean <- interaction.total/t.max
  return(interaction.mean)
}

### inputs the coordinates, edge lists and command to show or deny network
### returns a graph object used for metric calculations
get.network <- function(coords, edge.list, show.network){
  g <- graph.edgelist(edge.list[,1:2], directed=FALSE)                        # include g as weigted edgelist
  g <- set.edge.attribute(g, "weight", index=E(g), value=edge.list[,3])       # include weight
  w <- get.edge.attribute(g, "weight")                                        # make edge attribute 
  coords <- cbind(as.numeric(coords[,2]), as.numeric(coords[,1]))             # reverse x and y here
  modules <- leading.eigenvector.community(g)                                 # split communities using eigenvector algorithm
  if(show.network == TRUE){                                                   
    k <- tkplot(g, edge.width=w, vertex.color=membership(modules))            # use tkplot to visualise network, and show communities by colour
    tk_set_coords(k, coords=coords)                                           # place each individual using relative coordinates 
    tk_rotate(k, degree=-90)                                                  # rotate 90 degrees left for same view as spatial field
  }
  return(list("g"=g, "modules"=modules))
}

### function to report the progress of any simulation as it runs directly to the console
get.progress <- function(time.5, time.sim, time.start){
  time.elapsed <- round(proc.time()[3]-time.start, 1)
  time.5 <- binhf::shift(time.5, 1)
  time.5[1] <- proc.time()[3]-time.sim
  time.remaining <- round((n.sim-sim)*mean(time.5, na.rm=TRUE), 0)
  writeLines(paste("Completed ", sim, " of ", n.sim, " simulations in ", seconds_to_period(time.elapsed), ". ",
                   "Estimated time remaining: ", seconds_to_period(time.remaining), ".", sep=""))
  return(time.5)
}

### a function to calculate the confidence intervals of a regression model
ci.lines <- function(model){
  xm <- sapply(model[[12]][2], mean)
  n <- sapply(model[[12]][2], length)
  ssx <- sum(model[[12]][2]^2)-sum(model[[12]][2])^2/n
  s.t <- qt(0.975, (n-2))
  xv <- seq(min(model[[12]][2]), max(model[[12]][2]), length=100)
  yv <- coef(model)[1]+coef(model)[2]*xv
  se <- sqrt(summary(model)[[6]]^2*(1/n+(xv-xm)^2/ssx))
  ci <- s.t*se
  uyv <- yv+ci
  lyv <- yv-ci
  lines(xv, uyv, lty=2, col="black")
  lines(xv, lyv, lty=2, col="black")
}

## a function to calculate prediction intervals of a regression model
pr.lines <- function(metric, model){
  pr.int <- predict(model,interval="prediction")
  pr.lw <- pr.int[,2]
  pr.up <- pr.int[,3]
  # points(metric, pr.up, pch=20)
  # points(metric, pr.lw, pch=20)
  delta.y <- pr.up[1]-pr.up[2]
  delta.x <- metric[1]-metric[2]
  m <- delta.y/delta.x
  if (m > 0){
    decreasing <- FALSE
  } else{
    decreasing <- TRUE
  }
  lines(sort(metric), sort(pr.up, decreasing=decreasing))
  lines(sort(metric), sort(pr.lw, decreasing=decreasing))
}

## please comment out ci.lines or pi.lines if and when necessary
plot.metrics <- function(d, metric.names, metrics, y){
  if (y == "Lag"){y.number <- length(metric.names)+1}                     # which column is required (lag, post lag or rate)?
  else if (y == "Post Lag"){y.number <- length(metric.names)+2}
  else if (y == "Rate"){y.number <- length(metric.names)+3} 
  par(mfrow=c(3,5))                                                       # display graph 3 x 5 in plotting window
  for(i in 1:length(metric.names)){
    plot(metrics[, i, toString(d)], metrics[, y.number, toString(d)],     # plot metrics by y
         xlab=metric.names[i], 
         ylab=y,
         main = d)
    model <- lm(metrics[, y.number, toString(d)]~metrics[, i, toString(d)]) # choose a regression model to fit expected metric graphs
    abline(model)                                                         # Line of best fit, based on chosen regression model
    ci.lines(model)                                                       # Confidence interval
    # pr.lines(metrics[,i], model)                                        # Prediction interval
  }
}

### calculate the infection step
get.infection.step <- function(infection.status, infection.step, n, t){
  for (i in 1:n){
    if (infection.step[i] == 0 & infection.status[i] == 1){
      infection.step[i] <- t
    }
  }
  return(infection.step)
}

### creates a folder in working directory to save metrics
export.results <- function(metrics, metrics.leper){
  dir.create("Results")
  folder <- paste("Results/", cluster,"_","leper","_", n, "_", n.sim, sep="")
  dir.create(folder)
  write.table(metrics, paste(folder, "/metrics.csv", sep=""), row.names=FALSE, sep=",")
  write.table(metrics.leper,  paste(folder, "/metrics_leper.csv", sep=""), row.names=FALSE, sep=",")
}


### ------------------------------------------------------------------------------------------------------------ ###
###                                                 SET UP                                                       ###
### ------------------------------------------------------------------------------------------------------------ ###

metric.names <- c("Flow Betweenness", "Eigenvector", "Infocent", "Efficiency",                            # list of metric names (for plotting purposes)
                  "Modularity", "Strength", "Diameter", "Closeness",  "Betweenness", "Edge Betweenness", 
                  "Shortest Paths", "Weighted Betweenness", "Spreading Centrality", "Reach")
y.names <- c("Lag", "Post Lag", "Rate")                                                                    # list of infection spread measurements
metrics <- set.metrics(dispersion, metric.names, n.sim, y.names)                                           # 3d matrix for results
metrics.leper <- set.metrics(dispersion, metric.names, n.sim, y.names)
par(mfrow=c(1,2))                                                                                          # set graph format to 1 x 2
time.5 <- replace(vector(mode="numeric", length=0), 1:5, NA)                                               # set vector to store time taken for last 5 simulations
time.start <- proc.time()[3]                                                                               # record start time

### ------------------------------------------------------------------------------------------------------------ ###
###                                                   RUN                                                        ###
### ------------------------------------------------------------------------------------------------------------ ###

for(d in dispersion){                                            # d is any one dispersion value, 
  
  for(sim in 1:n.sim){
    time.sim <- proc.time()[3]                                   # record start time of this simulation
    coords <- get.coords(cluster, d, n, n.parents)               # set spatial field
    toroidal.distance <- get.toroidal.distance(coords)           # get Euclidian distance matrix
    
    ### RESET SOME MATRICES AND VECTORS FOR NEW SIMULATION ###
    interaction.total <- matrix(0, nrow=n, ncol=n)               # set a matrix to track total interactions
    infection.mean <- vector(mode="numeric", length=0)           # set a vector to track the mean infection status
    
    ### PICK THE LEPER ### 
    infection.status <- numeric(n)                               # set the infection status of all badgers to zero
    leper <- sample(1:n, 1)                                      # select a leper at random
    infection.status[leper] <- 1                                 # set the infection status of the leper to 1
    infection <- mean(infection.status)                          # calculate the initial mean infection status (1/n)
    infection.step <- numeric(n)                                 #
    t <- 0
    while(infection < infection.constant){
      t <- t+1                                                                               #
      interaction.expected <- get.interaction.expected(toroidal.distance, n)                 # get the expected interaction matrix for this time step     
      interaction.total <- get.interaction.total(interaction.expected, interaction.total)    # track the total number of interactions 
      infection.status <- get.infection.status(interaction.expected, n, p)                   # determine the new infected status 
      infection.step <- get.infection.step(infection.status, infection.step, n, t)           #
      infection <- mean(infection.status)                                                    # infection = mean(infected status)
      infection.mean <- c(infection.mean, infection)                                         # store infection
    }
    
    t.max <- length(infection.mean)                                              # determine the final time step
    interaction.mean <- get.interaction.mean(interaction.total, t.max)           # get the mean interaction matrix
    edge.list <- get.edge.list(interaction.mean, n)                              # get a weighted edgelist using mean interactions
    
    metrics[sim, "Lag", toString(d)] <- get.lag(infection.mean)                  # store lag
    metrics[sim, "Post Lag", toString(d)] <- get.postlag(infection.mean)         # store post lag
    metrics[sim, "Rate", toString(d)] <- get.rate(t.max)                         # store rate
    metrics.leper[sim, "Lag", toString(d)] <- get.lag(infection.mean)
    metrics.leper[sim, "Post Lag", toString(d)] <- get.postlag(infection.mean)
    metrics.leper[sim, "Rate", toString(d)] <- get.rate(t.max)
    
    if(show.spatial.field == TRUE){                                              # plot the spatial field and mean infection status with time
      plot(coords[,1], coords[,2], xlim=c(0,1), ylim=c(0,1), 
           main="Spatial Field",
           xlab="x-coordinate", 
           ylab="y-coordinate")
      plot(c(1:t.max), infection.mean, pch=20,
           main="Probability of Infection",
           xlab="Time", 
           ylab="Probability of Infection")
    }
    
    graph <- get.network(coords, edge.list, show.network)                        # get the network as a graph
    
    ### ------------------------------------------------------------------------------------------------------------ ###
    ###                                            GET AND STORE METRICS                                             ###
    ### ------------------------------------------------------------------------------------------------------------ ###
    
    ### SNA METRICS ###
    
    ### note that metrics labelled as population-level will not be calculated for leper individual statisitics
    
    
    fb <- flowbet(interaction.mean)
    metrics[sim, "Flow Betweenness", toString(d)] <- mean(fb)
    metrics.leper[sim, "Flow Betweenness", toString(d)] <- fb[leper]
    
    ev <- evcent(interaction.mean)
    metrics[sim, "Eigenvector", toString(d)] <- mean(ev)
    metrics.leper[sim, "Eigenvector", toString(d)] <- ev[leper]
    
    ic <- infocent(interaction.mean)
    metrics[sim, "Infocent", toString(d)] <- mean(ic)
    metrics.leper[sim, "Infocent", toString(d)] <- ic[leper]
    
    metrics[sim, "Efficiency", toString(d)] <- efficiency(interaction.mean)
    metrics.leper[sim, "Efficiency", toString(d)] <- NA                            # efficiency is a population level metric
    
    ### IGRAPH METRICS ###
    metrics[sim, "Modularity", toString(d)] <- modularity(graph$modules)
    metrics.leper[sim, "Modularity", toString(d)] <- NA                            # modularity is a population level metric
    
    sg <- strength(graph$g)
    metrics[sim, "Strength", toString(d)] <- mean(sg)
    metrics.leper[sim, "Strength", toString(d)] <- sg[leper]
    
    metrics[sim, "Diameter", toString(d)] <- diameter(graph$g)
    metrics.leper[sim, "Diameter", toString(d)] <- NA                              # diameter is a population level metric
    
    cs <- closeness.estimate(graph$g, cutoff = 5)
    metrics[sim, "Closeness", toString(d)] <- mean(cs)
    metrics.leper[sim, "Closeness", toString(d)] <- cs[leper]
    
    bw <- betweenness.estimate(graph$g, cutoff = 5)
    metrics[sim, "Betweenness", toString(d)] <- mean(bw)
    metrics.leper[sim, "Betweenness", toString(d)] <- bw[leper]
    
    eb <- edge.betweenness(graph$g)
    metrics[sim, "Edge Betweenness", toString(d)] <- mean(eb)
    metrics.leper[sim, "Edge Betweenness", toString(d)] <- NA
    
    sp <- shortest.paths(graph$g)
    metrics[sim, "Shortest Paths", toString(d)] <- mean(sp)
    metrics.leper[sim, "Shortest Paths", toString(d)] <- mean(sp[leper,])
    
    wb <- betweenness_w(interaction.mean, directed=NULL, alpha=1)[,2]
    metrics[sim, "Weighted Betweenness", toString(d)] <- mean(wb)
    metrics.leper[sim, "Weighted Betweenness", toString(d)] <- wb[leper]
    
    #### CePa ###
    spd <- spread(graph$g)
    metrics[sim, "Spreading Centrality", toString(d)] <- mean(spd)
    metrics.leper[sim, "Spreading Centrality", toString(d)] <- spd[leper]
    
    rh <- reach(graph$g, weights=E(graph$g)$weight)
    metrics[sim, "Reach", toString(d)] <- mean(rh)
    metrics.leper[sim, "Reach", toString(d)] <- rh[leper]
    
    time.5 <- get.progress(time.5, time.sim, time.start)                              # Print progress in console
  }
}

  ### ------------------------------------------------------------------------------------------------------------ ###
  ###                                                 PLOT METRICS                                                 ###
  ### ------------------------------------------------------------------------------------------------------------ ###
  
if (show.metrics == TRUE){
  for (d in dispersion){
    for (y in y.names){
      plot.metrics(d, metric.names, metrics, y)
    }
  }
}

### ------------------------------------------------------------------------------------------------------------ ###
###                                                 EXPORT RESULTS                                               ###
### ------------------------------------------------------------------------------------------------------------ ###

if (save.results == TRUE){                            # save as csv file in working directory
  export.results(metrics, metrics.leper)
}
