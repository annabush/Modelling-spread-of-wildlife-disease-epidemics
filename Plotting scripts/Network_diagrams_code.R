# Network Plotter

# NB: Comment out/include appropriate lines specific to a leper or a probability network 

# network <- "observed"
network <- "expected"

par(mfrow=c(1,1))
weight.lim <- c(1.5,2.5) 
weight <- c(0.75, 1.5)
color <- c("lightgrey", "black")


range01 <- function(x){
  (x-min(x))/(max(x)-min(x))
  }


get.network <- function(adj){
  g <- graph.adjacency(adj, mode="undirected", weighted=TRUE)
  # SET EDGE WIDTH
  E(g)$width <- ifelse(E(g)$weight > weight.lim[2] , weight[2], ifelse(E(g)$weight > weight.lim[1], weight[1], 0))
  # SET EDGE COLOUR
  E(g)$color <- ifelse(E(g)$weight > weight.lim[2] , color[2], ifelse(E(g)$weight > weight.lim[1], color[1], 0))
  # SET VERTEX SHAPE
  V(g)$shape <- "circle"
  if (network == "observed")
  {
    V(g)$shape[leper] <- "square"                                     
  }
  # SET VERTEX COLOUR
  if (network == "observed")
  {
    V(g)$color <- ifelse(infection.status == 1, "tomato2", "royalblue2")     
  } else if (network = "expected")
  {
    V(g)$color <- rgb(range01(infection.status),0.36,1-range01(infection.status))
  }
  # # SET VERTEX SIZE
  V(g)$size <- 3
  # LABEL 
  if (network == "expected")
  {
    V(g)$label <- round(bw,1)
  } else if (network == "observed")
  {
    V(g)$label <- round(fb,1) 
  }
  # SET LABEL SIZE
  V(g)$label.cex <- 2
  return(g)
}


g <- get.network(interaction.mean)

# Plot graph
plot(g, 
     layout=as.matrix(coords), 
     vertex.label.color="black", 
     vertex.label.dist=0.4)

# Set font
op <- par(family = "serif")

# LEGENDs
if (network == "observed")
{
  legend(0.3, 1,
         legend=c("Infected", "Uninfected", "Leper"),
         col=c("tomato2", "royalblue2", "tomato2"),
         y.intersp=0.3,
         x.intersp=0.5,
         pch=c(19,19,15), 
         pt.cex=4,
         cex=2,
         bty="n")
} else if (network == "expected")
{
  legend(0.3,1,
         legend=c(paste("Probability of infection = ", round(min(infection.status),2), sep=""), 
                  paste("Probability of infection = ", round(max(infection.status),2), sep="")),
         col=c("royalblue2", "tomato2"),
         y.intersp=0.3,
         x.intersp=0.5,
         pch=c(19,19), 
         pt.cex=4,
         cex=2,
         bty="n")
}

# Apply font
par(op)




