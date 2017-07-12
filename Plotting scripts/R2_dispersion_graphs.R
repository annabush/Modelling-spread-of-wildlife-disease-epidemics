disp.names <- c("0.1", "0.125", "0.15", "0.2", "0.3", "1")

# network <- "expected"
network <- "observed"

col <- c("royalblue2", "tomato2", "green3")
lwd <- 2
if (network == "observed")
{
  ylim <- c(0.45, 1)        # y limits for graph - OBSERVED
} else if (network == "expected")
{
  ylim <- c(0, 0.6)         # y limits for graph - EXPECTED
}
y.intersp <- 1.4          # y spacing for legend
x.space <- 2.5            # x spacing for legend
point <- 1.5              # point size
label <- 1.2              # axis label size
leg.lab <- 1              # legend label size
axis <- 1                 # axis size

line <- function(y, col){
  for (i in 1:5){
    segments(i, y[i], i+1, y[i+1], col=col, lwd=lwd)
  }
}

if (network == "observed")
{
  # OBSERVED
  obs <- matrix(nrow=6, ncol=3)
  colnames(obs) <- c("Bet_r", "Bet_l", "Bet_p")
  rownames(obs) <- disp.names
  obs <- as.data.frame(obs)
  obs$Bet_r <- c(0.3293, 0.2838, 0.2509, 0.3545, 0.4599, 0.4685)
  obs$Bet_l <- c(0.07812, 0.05462, 0.0691, 0.1269, 0.1101, 0.1544)
  obs$Bet_p <- c(0.2743, 0.2734, 0.1613, 0.275, 0.3145, 0.2656)
  # OBSERVED - LEPER
  lep <- matrix(nrow=6, ncol=3)
  colnames(lep) <- c("Str_l", "Bet_r", "Bet_p")
  rownames(lep) <- disp.names
  lep <- as.data.frame(lep)
  lep$Bet_r <- c(0.016, 0.0328, 0.06016, 0.009354, 0.009404, 0.007343)
  lep$Str_l <- c(0.02884, 0.1029, 0.06085, 0.0469, 0.0541, 0.01961)
  lep$Bet_p <- c(0.01192, 0.0009899, 0.02723, 0.01358, 0.0135, 0.008741)
} else if  (network == "expected")
{
  exp <- matrix(nrow=6, ncol=6)
  colnames(exp) <- c("Eff_r", "Eff_l", "Eff_p", "Flo_r", "Flo_l", "Flo_p")
  rownames(exp) <- disp.names
  exp <- as.data.frame(exp)
  exp$Eff_r <- c(0.9963, 0.9963, 0.9965, 0.9975, 0.9978, 0.9976)
  exp$Eff_l <- c(0.9952, 0.9951, 0.9943, 0.9962, 0.9927, 0.976)
  exp$Eff_p <- c(0.9954, 0.9953, 0.9957, 0.9966, 0.9969, 0.9963)
  exp$Flo_r <- c(0.8167, 0.7974, 0.8215, 0.8844, 0.8515, 0.4733)
  exp$Flo_l <- c(0.8322, 0.8145, 0.8357, 0.8976, 0.8639, 0.4727)
  exp$Flo_p <- c(0.8124, 0.7927, 0.8173, 0.8805, 0.8471, 0.4706)
}



par(mar=c(5,5,1,1)) 
# Draw empty plot
plot(1, type="n",  xaxt="n",
     xlab="Dispersion", 
     ylab=expression(R^2), 
     xlim=c(1,6), 
     ylim=ylim, 
     cex.lab=label,
     cex.axis=axis)

# Change x-axis
axis(side = 1,
     at = 1:6,
     labels = disp.names,
     cex.axis = axis)


if (network == "expected")
{
  # ADD POINTS
  points(1:6, exp$Eff_r, col=col[1], pch=1, cex=point)
  points(1:6, exp$Eff_l, col=col[2], pch=1, cex=point)
  points(1:6, exp$Eff_p, col=col[3], pch=1, cex=point)
  points(1:6, exp$Flo_r, col=col[1], pch=2, cex=point)
  points(1:6, exp$Flo_l, col=col[2], pch=2, cex=point)
  points(1:6, exp$Flo_p, col=col[3], pch=2, cex=point)
  # ADD LINES
  line(exp$Eff_r, col[1])
  line(exp$Eff_l, col[2])
  line(exp$Eff_p, col[3])
  line(exp$Flo_r, col[1])
  line(exp$Flo_l, col[2])
  line(exp$Flo_p, col[3])
  # LEGEND
  legend(x.space, 0.6, 
       c("Efficiency", "Flow-betweenness"),
       bty="n",
       pch=c(1,2),
       cex=leg.lab,
       pt.cex = point,
       y.intersp=y.intersp)
} else if (network == "observed")
{
  # ADD POINTS
  points(1:6, obs$Bet_r, col=col[1], pch=1, cex=1.2)
  points(1:6, obs$Bet_l, col=col[2], pch=1, cex=1.2)
  points(1:6, obs$Bet_p, col=col[3], pch=1, cex=1.2)
  points(1:6, lep$Bet_r, col=col[1], pch=3, cex=1.2)
  points(1:6, lep$Str_l, col=col[2], pch=2, cex=1.2)
  points(1:6, lep$Bet_p, col=col[3], pch=3, cex=1.2)
  # ADD LINES
  line(obs$Bet_r, col[1])
  line(obs$Bet_l, col[2])
  line(obs$Bet_p, col[3])
  line(lep$Bet_r, col[1])
  line(lep$Str_l, col[2])
  line(lep$Bet_p, col[3])
  # LEGEND
  legend(x.space, 0.62, 
         c("Betweenness (individual-based)", "Strength (source-infector analyses)", "Betweenness (source-infector analyses)"),
         bty="n",
         pch=c(1,2,3),
         pt.cex=point,
         y.intersp=y.intersp)
}


# LEGEND
legend(1, 0.62, 
       c("Rate", "Lag", "Post-lag"),
       bty="n",
       lty=c(1,1,1),
       lwd=c(lwd,lwd,lwd),
       cex=leg.lab,
       col=col,
       y.intersp=y.intersp)

