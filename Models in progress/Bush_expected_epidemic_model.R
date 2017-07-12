### ------------------------------------------------------------------------------------------------------------ ###
### Author: Anna F. Bush                                                                                         ###
### Date last revised: 15th July, 2016                                                                           ###
### ------------------------------------------------------------------------------------------------------------ ###
### 
###                              EXPECTED MODEL SHOWING PROBABALISTIC DISEASE DYNAMICS
###
### ------------------------------------------------------------------------------------------------------------ ###
### A program to simulate disease dynamics within a network, based on an increasing probability 
### of disease transmission as infection persists.
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
dispersion <- c(0.1)
cluster <- "uniform"                            ### normal or uniform                                 
n <- 20                                         ### number of individuals                        
p <- 1*10^-3                                    ### infectiveness of disease                     
infection.constant <- 0.5                       ### percentage of infected individuals at cut off   
n.sim <- 50                                     ### number of model iterations         
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
  metrics <- array(NA,                                                                  # NA to return error if value not filled 
                   dim=c(n.sim, length(metric.names)+3, length(dispersion)),            # n.sim x metric.names x dispersion list
                   dimnames=list(1:n.sim, c(metric.names, y.names), dispersion))        # Name the dimensions
  return(metrics)                                                                       # Return as data a data frame
}

### generates a random spatial field using rpcp function
### given the cluster type, dispersion value and number of individuals
### cartesian coordinates are returned
get.coords <- function(cluster, dispersion, n, n.parents){
  coords <- rpcp(nparents=n.parents, npoints=n, lambda=0.5, nsim=1, cluster=cluster, dispersion=dispersion, infection=FALSE)
  coords <- cbind(as.numeric(coords$xyt[,1]), as.numeric(coords$xyt[,2]))
  colnames(coords)[1:2] <- c("x", "y")
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
      delta.x <- min(delta.x,1-delta.x)                # select the shortest route in the x direction
      delta.y <- min(delta.y,1-delta.y)                # select the shortest route in the y direction
      dist <- sqrt(delta.x^2+delta.y^2)                # Pythagoras
      toroidal.distance[i, j] <- dist                  # fill matrix
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

### keeps track of the ongoing score of interactions between each i and j over multiple iterations
get.interaction.total <- function(interaction.expected, interaction.total){
  interaction.total <- interaction.total + interaction.expected
  return(interaction.total)
}


### gives the probability that each individual will infect any other individual in network
### based on the chain binomial epidemic model equation 
get.infection.matrix <- function(infection.status, interaction.expected, n, p){
  infection.matrix <- matrix(0, nrow=n, ncol=n)
  for (i in 1:n){                                                   # for each row
    for (j in 1:n){                                                 # for each column
     infection.matrix[i, j] <- 1-(1-infection.status[j]*p)^interaction.expected[i,j] # p = 1-(1-p*sj)^ei,j
    }
  }
  return(infection.matrix)
}


### Calculates the probability that each individual has been infected in any one time step
get.infection.probability <- function(infection.matrix, n){
  infection.probability <- numeric(n)
  for (i in 1:n){                                  # for each row 
    p <- 1
    for (j in 1:n){                                # for each column
      p <- p*(1-infection.matrix[i, j])            # multiply all (1-infection probabilities) for this row
    }
    infection.probability[i] <- 1-p                # probability that each individual becomes infected in this timestep
  }
  return(infection.probability)
}


### generates a vector of probabilities that each individual is infected in each iteration of the simulation
get.infection.status <- function(infection.probability, infection.status){
  infection.status <- 1-(1-infection.status)*(1-infection.probability)
  return(infection.status)
}


### calculates lag
get.lag <- function(infection.mean, n){
  t2 <- min(which(infection.mean>2/n))-1                     # t1 = time steps to get above 2/n
  t1 <- t2-1                                                 # t2 = time step before infection gets to 2/n
  y1 <- infection.mean[t1+1]                                   # infection at t1
  y2 <- infection.mean[t2+1]                                   # infection at t2
  time <- t1 + ((2/n-y1)*(t2-t1)/(y2-y1))                    # linear interpolation 
  lag <- 1/(n*time)                                          # gradient = (2/n - 1/n)/(time - 0) = 1/(n*time)
  return(lag)
}


### calculates post lag
get.post.lag <- function(infection.constant, infection.mean, n, t.max){
  t2 <- min(which(infection.mean>2/n))-1                     # t1 = time steps to get above 2/n
  t1 <- t2-1                                                 # t2 = time step before infection gets to 2/n
  y1 <- infection.mean[t1+1]  
  y2 <- infection.mean[t2+1]
  time.1 <- t1 + ((2/n-y1)*(t2-t1)/(y2-y1))
  t1 <- t.max-1
  t2 <- t.max
  y1 <- infection.mean[t1+1]
  y2 <- infection.mean[t2+1]
  time.2 <- t1 + ((infection.constant-y1)/(y2-y1))
  post.lag <- (infection.constant-2/n)/(time.2-time.1)
  return(post.lag)
}


### calculates rate 
get.rate <- function(infection.constant, infection.mean, t.max){
  t1 <- t.max-1
  t2 <- t.max
  y1 <- infection.mean[t1+1]
  y2 <- infection.mean[t2+1]
  time <- t1 + ((infection.constant-y1)/(y2-y1))            # linear interpolation
  rate <- (infection.constant-1/n)/time
  return(rate)
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

### interaction mean is the total number of interactions 
### divided by the time in which these interactions took place
get.interaction.mean <- function(interaction.total, t.max){
  interaction.mean <- interaction.total/t.max
  return(interaction.mean)
}

### inputs the coordinates, edge lists and command to show or deny network
### returns a graph object used for metric calculations
get.network <- function(coords, edge.list, show.network){
  g <- graph.edgelist(edge.list[,1:2], directed=FALSE)                   # include g as weigted edgelist
  g <- set.edge.attribute(g, "weight", index=E(g), value=edge.list[,3])  # include weight
  w <- get.edge.attribute(g, "weight")                                   # make edge attribute
  coords <- cbind(as.numeric(coords[,2]), as.numeric(coords[,1]))        # reverse x and y
  modules <- leading.eigenvector.community(g)                            # split communities using eigenvector algorithm
  if(show.network == TRUE){                                          
    k <- tkplot(g, edge.width=w, vertex.color=membership(modules))       # can use tkplot to visualise network
    tk_set_coords(k, coords=coords)                                      # use relative coordinates
    tk_rotate(k, degree=-90)                                             # rotate 90 degrees left for same view as spatial field
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

### a function to calculate the prediction intervals of a regression model
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
  if (y == "Lag"){y.number = 1+length(metric.names)}                     # which column is required (lag, post lag or rate)?
  else if (y == "Post Lag"){y.number=2+length(metric.names)}
  else if (y == "Rate"){y.number = 3+length(metric.names)} 
  par(mfrow=c(3,5))                                                      # display graph 3 x 5 to include call metrics
  for(i in 1:length(metric.names)){
    plot(metrics[, i, toString(d)], metrics[, y.number, toString(d)], 
         xlab=metric.names[i], 
         ylab=y,
         main=d)
    model <- lm(metrics[, y.number, toString(d)]~metrics[, i, toString(d)])
    abline(model)
    ci.lines(model)                                                       # Confidence interval
    # pr.lines(metrics[,i], model)                                        # Prediction interval
  }
}


export.results <- function(metrics, cluster, n, n.sim){
  dir.create("Results")
  folder <- paste("Results/", cluster, "_","probability","_", n, "_", n.sim, sep="")
  dir.create(folder)
  write.table(metrics, paste(folder, "/metrics.csv", sep=""), row.names=FALSE, sep=",")
  }


### ------------------------------------------------------------------------------------------------------------ ###
###                                                 SET UP                                                       ###
### ------------------------------------------------------------------------------------------------------------ ###

metric.names <- c("Flow Betweenness", "Eigenvector", "Infocent", "Efficiency", 
                  "Modularity", "Strength", "Diameter", "Closeness",  "Betweenness", "Edge Betweenness", 
                  "Shortest Paths", "Weighted Betweenness", "Spreading Centrality", "Reach")
y.names <- c("Lag", "Post Lag", "Rate")
metrics <- set.metrics(dispersion, metric.names, n.sim, y.names)
par(mfrow=c(1,2))
time.5 <- replace(vector(mode="numeric", length=0), 1:5, NA)   # make a vector to store time values
time.start <- proc.time()[3]                                   # Record absolute start time

### ------------------------------------------------------------------------------------------------------------ ###
###                                                   RUN                                                        ###
### ------------------------------------------------------------------------------------------------------------ ###

for (d in dispersion){
  
  for(sim in 1:n.sim){
    time.sim <- proc.time()[3]                                    # Record start time of this simulation
    coords <- get.coords(cluster, d, n, n.parents)                # Generate coordinates of individuals
    toroidal.distance <- get.toroidal.distance(coords)            # Calculate the Euclidian distance matrix
    
    ### RESET SOME MATRICES AND VECTORS FOR NEW SIMULATION ###
    interaction.total <- matrix(0, nrow=n, ncol=n)                # Tracks the total numer of interactions for each individual
    infection.status <- rep(1/n, n)                               # Vector of probabilities that each individual is infected
    infection.mean <- 1/n
  
    while(mean(infection.status) < infection.constant){
      interaction.expected <- get.interaction.expected(toroidal.distance, n)                            # Generate the expected interactions between individuals
      interaction.total <- get.interaction.total(interaction.expected, interaction.total)               # Update the total number of interactions
      
      infection.matrix <- get.infection.matrix(infection.status, interaction.expected, n, p)            # Return matrix if infection probabilities
      infection.probability <- get.infection.probability(infection.matrix, n)                           # Calculate the probability that each individual has been infected in this timestep
      infection.status <- get.infection.status(infection.probability, infection.status)                 # Calculate the probability that each individual is now infected
      
      infection.mean <- append(infection.mean, mean(infection.status))                                  # Store the mean infection.status
    }
    
    t.max <- length(infection.mean)-1                                              # t starts at 0 so do length -1
    
    if(show.spatial.field == TRUE){                                                # Plot coordinates and infection.status
      plot(coords[,1], coords[,2], xlim=c(0,1), ylim=c(0,1), 
           main="Spatial Field",
           xlab="x-coordinate", 
           ylab="y-coordinate")
      plot(c(0:t.max), infection.mean, pch=20,
           main="Probability of Infection",
           xlab="Time",
           ylab="Probability of Infection")
    }
    
    interaction.mean <- get.interaction.mean(interaction.total, t.max)    
    edge.list <- get.edge.list(interaction.mean, n)
    graph <- get.network(coords, edge.list, show.network)
    
    metrics[sim, "Lag", toString(d)] <- get.lag(infection.mean, n)
    metrics[sim, "Post Lag", toString(d)] <- get.post.lag(infection.constant, infection.mean, n, t.max)
    metrics[sim, "Rate", toString(d)] <- get.rate(infection.constant, infection.mean, t.max)
    
    ### ------------------------------------------------------------------------------------------------------------ ###
    ###                                            GET AND STORE METRICS                                             ###
    ### ------------------------------------------------------------------------------------------------------------ ###
    
    ### SNA METRICS ###
    metrics[sim, "Flow Betweenness", toString(d)] <- mean(flowbet(interaction.mean))
    metrics[sim, "Eigenvector", toString(d)] <- mean(evcent(interaction.mean))
    metrics[sim, "Infocent", toString(d)] <- mean(infocent(interaction.mean))
    metrics[sim, "Efficiency", toString(d)] <- efficiency(interaction.mean)
    ### IGRAPH METRICS ###
    metrics[sim, "Modularity", toString(d)] <- modularity(graph$modules)
    metrics[sim, "Strength", toString(d)] <- mean(strength(graph$g))
    metrics[sim, "Diameter", toString(d)] <- diameter(graph$g)
    metrics[sim, "Closeness", toString(d)] <- mean(closeness.estimate(graph$g, cutoff = 5))
    metrics[sim, "Betweenness", toString(d)] <- mean(betweenness.estimate(graph$g, cutoff = 5))
    metrics[sim, "Edge Betweenness", toString(d)] <- mean(edge.betweenness(graph$g))
    metrics[sim, "Shortest Paths", toString(d)] <- mean(shortest.paths(graph$g))
    metrics[sim, "Weighted Betweenness", toString(d)] <- mean(betweenness_w(interaction.mean, directed=NULL, alpha=1)[,2])
    ####CePa
    metrics[sim, "Spreading Centrality", toString(d)] <- mean(spread(graph$g))
    metrics[sim, "Reach", toString(d)] <- mean(reach(graph$g, weights=E(graph$g)$weight))
    
    time.5 <- get.progress(time.5, time.sim, time.start)                              # Print progress
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

if (save.results == TRUE){
  export.results(metrics, cluster, n, n.sim)
}