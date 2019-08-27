inputFolder="../data/calibration_2019_small/CENT"
inputIr="$inputFolder/calibration_2019_00_C_20190509_035506.681605_ir.tif"
inputOptic="$inputFolder/debayer/calibration_2019_00_C_20190509_035506.681605_rgb.tif"

outputFolder="results_ir_to_optic_CENT_calibration_2019_00_C_20190509_035506.681605"
outputIrPhase="$outputFolder/ir_phase.mhd"
outputIrPhaseIn="$outputFolder/ir.mhd"
outputIrPoints="$outputFolder/ir_points"
outputOpticPhase="$outputFolder/optic_phase.mhd"
outputOpticPhaseIn="$outputFolder/optic.mhd"
outputOpticPoints="$outputFolder/optic_points"

mkdir $outputFolder


echo "Phase"
./0_PhaseSymmetry $inputIr 1  $outputIrPhase $outputIrPhaseIn
./0_PhaseSymmetry $inputOptic 0  $outputOpticPhase $outputOpticPhaseIn

echo "Point Set"
./1_NarrowBandPointSet $outputIrPhase $outputIrPoints 0.8
./1_NarrowBandPointSet $outputOpticPhase $outputOpticPoints 8

#Trimmed Euclidean Registration
trimmedTransform="$outputFolder/trimmed-euclidean.h5"
./2_PointSetRegistration "$outputIrPoints.gii" "$outputOpticPoints.gii" $trimmedTransform 3
./3_TransformPointSetAndImage $trimmedTransform "$outputIrPoints.gii" "$outputIrPoints-trimmed-euclidean.off" \
                             "$outputOpticPoints.gii" "$outputOpticPoints-trimmed-euclidean.off" \
                              $inputIr  $inputOptic $outputIrPhase $outputOpticPhase "$outputFolder/trimmed-euclidean-viame.h5"

#Euclidean Registration
euclideanTransform="$outputFolder/euclidean.h5"
#./2_PointSetRegistration "$outputIrPoints.gii" "$outputOpticPoints.gii" $euclideanTransform 0
#./3_TransformPointSetAndImage $euclideanTransform "$outputIrPoints.gii" "$outputIrPoints-euclidean.off" \
#                              "$outputOpticPoints.gii" "$outputOpticPoints-euclidean.off" \
#                              $inputIr  $inputOptic "$outputFolder/euclidean-viame.h5"

#Optimal Transport Registration
otTransform="$outputFolder/ot.h5"
#./2_PointSetRegistration "$outputIrPoints.gii" "$outputOpticPoints.gii" $otTransform 5
#./3_TransformPointSetAndImage $otTransform "$outputIrPoints.gii" "$outputIrPoints-ot.off" \
#                              "$outputOpticPoints.gii" "$outputOpticPoints-ot.off" \
#                              $inputIr  $inputOptic "$outputFolder/ot-viame.h5"


