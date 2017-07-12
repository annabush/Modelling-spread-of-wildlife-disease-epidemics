# CODE TO PLOT METRICS AGAINST RATE LAG AND POST LAG


# SELECT FILE PATH TO RESULTS
# file.path <- "Documents/Anna/Research Project/Results/Probability_metrics_uniform_250.csv"
file.path <- "Documents/Anna/Research Project/Results/Leper_metrics_uniform_250.csv"
# file.path <- "Documents/Anna/Research Project/Results/Leper_individual_uniform_250.csv"

# SELECT graphs per row and column 
nrow <- 1
ncol <- 1

# SELECT DISPERSIONS
# disp <- 0.1                          # Plot all graphs for chosen dispersions
# disp <- unique(results$Dispersion)   # Plots graphs for all dispersions
disp <- "All Dispersions"              # All dispersions on one graph

# SELECT LABEL SIZE
label <- 1.2

results <- read.csv(file.path)

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

for (d in disp){
  if (is.numeric(d)){
    df <- results[results$Dispersion == d,]
  } else {df <- results}
  
  x.var <- list("Flow Betweenness"=df$Flow.Betweenness, "Eigenvector"=df$Eigenvector, "Information centrality"=df$Infocent, 
                "Efficiency"=df$Efficiency, "Modularity"=df$Modularity, "Strength"=df$Strength, "Diameter"=df$Diameter, 
                "Closeness"=df$Closeness, "Betweenness centrality"=df$Betweenness, "Edge Betweenness"=df$Edge.Betweenness, 
                "Shortest Paths"=df$Shortest.Paths, "Weighted Betweenness"=df$Weighted.Betweenness, 
                "Spreading Centrality"=df$Spreading.Centrality, "Reach"=df$Reach)
  y.var <- list("Rate"=df$Rate, "Lag"=df$Lag,"Post Lag"=df$Postlag)
  
  par(mfrow=c(nrow,ncol))
  for (j in 1:length(y.var)){
    for (i in 1:length(x.var)){
      if (is.na(x.var[[i]][1]) == FALSE){
        plot(x.var[[i]], y.var[[j]],
             xlab=names(x.var[i]),
             ylab=names(y.var[j]),
             main= paste("Dispersion = ",d,sep=""),
             cex.lab=label)
        model <- lm(y.var[[j]]~x.var[[i]])
        abline(model)
        ci.lines(model)
      }
    }
  }
}