outputFolder <- "./results_ir_to_optic_CENT_calibration_2019_00_C_20190509_035506.681605"
ir    <- read.table( sprintf("%s/ir_points.off", outputFolder), header=FALSE, skip=2)
optic <- read.table( sprintf("%s/optic_points.off", outputFolder), header=FALSE, skip=2)


plot( optic[, 1:2], pch=".")
points( ir[,1:2], col="#0000FF22", pch=".")


dev.new()
ir.trimmed <- read.table( sprintf("%s/ir_points-trimmed-euclidean.off", outputFolder), header=FALSE, skip=2)
plot( optic[, 1:2], pch=".")
points( ir.trimmed[,1:2], col="#FF00FF88", pch=".")
title("Trimmed Euclidean")

dev.new()
optic.trimmed <- read.table( sprintf("%s/optic_points-trimmed-euclidean.off", outputFolder), header=FALSE, skip=2)
plot( ir[, 1:2], pch=".")
points( optic.trimmed[,1:2], col="#FF00FF88", pch=".")
title("Trimmed Euclidean Inverse")


#dev.new()
#ir.euclidean <- read.table( sprintf("%s/ir_points-euclidean.off", outputFolder), header=FALSE, skip=2)
#plot( optic[, 1:2], pch=".")
#points( ir[,1:2], col="#0000FF22", pch=".")
#points( ir.euclidean[,1:2], col="#FF00FF88", pch=".")
#title("Euclidean")


#dev.new()
#ir.ot <- read.table( sprintf("%s/ir_points-ot.off", outputFolder), header=FALSE, skip=2)
#plot( optic[, 1:2], pch=".")
#points( ir[,1:2], col="#0000FF22", pch=".")
#points( ir.ot[,1:2], col="#FF00FF88", pch=".")
#title("OT")
