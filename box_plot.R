set.seed(1)
n <- 1200
dat <- data.frame(
  x = gl(n=4, k=n/4),
  y = rnorm(n)
)
myColours = c(1, "steelblue", "#FFBB00", rgb(0.4, 0.2, 0.3))
#myColoursAlpha <- add.alpha(myColours, alpha=0.4)
op <- par(mfrow=c(1,2), mar=c(2,2,3,1))
boxplot(y ~ x, data=dat, outline=FALSE,
        axes=FALSE, main="alpha=1")
points(x=jitter(as.numeric(dat$x)), y=dat$y, 
       col=myColours[dat$x], pch=19)
