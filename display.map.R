ir <- read.csv("fixedTxfTrimmed.csv", header=F)
optic <- read.csv("movingTrimmed.csv", header=F)
map <-read.csv("tmp.csv", header=F)

plot(optic, pch=".")
cols = rgb(0,0,0, alpha=map[,3] / max(map[,3])*0.1)
segments( x0=ir[map[,1]+1,1], y0=ir[map[,1]+1,2], x1=optic[map[,2]+1,1], y1=optic[map[,2]+1,2], col=cols)
points(ir, col="red", pch=".")

