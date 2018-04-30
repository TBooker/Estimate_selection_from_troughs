../../../poly_dfe/polyDFE/polyDFE -d analysis/Ns400_pa_0.000125_polyDFE.txt -r ../../../poly_dfe/polyDFE/input/init_model_BandC.txt 1 -r ../../../poly_dfe/polyDFE/input/range_model_BandC.txt 1000 -m B > analysis/Ns400_pa_0.000125_polyDFE.fullDFE.withDiv.out &
../../../poly_dfe/polyDFE/polyDFE -d analysis/Ns400_pa_0.000125_polyDFE.txt -r ../../../poly_dfe/polyDFE/input/init_model_BandC.txt 2 -r ../../../poly_dfe/polyDFE/input/range_model_BandC.txt 1000 -m B > analysis/Ns400_pa_0.000125_polyDFE.dDFE.withDiv.out &

../../../poly_dfe/polyDFE/polyDFE -d analysis/Ns400_pa_0.000125_polyDFE.txt -r ../../../poly_dfe/polyDFE/input/init_model_BandC.txt 1 -r ../../../poly_dfe/polyDFE/input/range_model_BandC.txt 1000 -m B -w > analysis/Ns400_pa_0.000125_polyDFE.fullDFE.noDiv.out &
../../../poly_dfe/polyDFE/polyDFE -d analysis/Ns400_pa_0.000125_polyDFE.txt -r ../../../poly_dfe/polyDFE/input/init_model_BandC.txt 2 -r ../../../poly_dfe/polyDFE/input/range_model_BandC.txt 1000 -m B -w > analysis/Ns400_pa_0.000125_polyDFE.dDFE.noDiv.out &

../../../poly_dfe/polyDFE/polyDFE -d analysis/Ns400_pa_0.000125_polyDFE.txt -r ../../../poly_dfe/polyDFE/input/range_model_BandC.txt 1000 -m B -w > analysis/Ns400_pa_0.000125_polyDFE.BIGRANGE.noDiv.out &

../../../poly_dfe/polyDFE/polyDFE -d analysis/Ns200_pa_0.00025_polyDFE.txt -r ../../../poly_dfe/polyDFE/input/range_model_BandC.txt 1000 -m B > analysis/Ns200_pa_0.00025_polyDFE.BIGRANGE.withDiv.out &
../../../poly_dfe/polyDFE/polyDFE -d analysis/Ns200_pa_0.00025_polyDFE.txt  -r ../../../poly_dfe/polyDFE/input/range_model_BandC.txt 1000 -m B -w > analysis/Ns200_pa_0.00025_polyDFE.BIGRANGE.noDiv.out &

../../../poly_dfe/polyDFE/polyDFE -d analysis/Ns100_pa_0.0005_polyDFE.txt -r ../../../poly_dfe/polyDFE/input/range_model_BandC.txt 1000 -m B > analysis/Ns100_pa_0.0005_polyDFE.BIGRANGE.withDiv.out &
../../../poly_dfe/polyDFE/polyDFE -d analysis/Ns100_pa_0.0005_polyDFE.txt -r ../../../poly_dfe/polyDFE/input/range_model_BandC.txt 1000 -m B -w > analysis/Ns100_pa_0.0005_polyDFE.BIGRANGE.noDiv.out &

../../../poly_dfe/polyDFE/polyDFE -d analysis/Ns50_pa_0.001_polyDFE.txt -r ../../../poly_dfe/polyDFE/input/range_model_BandC.txt 1000 -m B > analysis/Ns50_pa_0.001_polyDFE.BIGRANGE.withDiv.out &
../../../poly_dfe/polyDFE/polyDFE -d analysis/Ns50_pa_0.001_polyDFE.txt -r ../../../poly_dfe/polyDFE/input/range_model_BandC.txt 1000 -m B -w > analysis/Ns50_pa_0.001_polyDFE.BIGRANGE.noDiv.out &

../../../poly_dfe/polyDFE/polyDFE -d analysis/Ns25_pa_0.002_polyDFE.txt -r ../../../poly_dfe/polyDFE/input/range_model_BandC.txt 1000 -m B > analysis/Ns25_pa_0.002_polyDFE.BIGRANGE.withDiv.out &
../../../poly_dfe/polyDFE/polyDFE -d analysis/Ns25_pa_0.002_polyDFE.txt -r ../../../poly_dfe/polyDFE/input/range_model_BandC.txt 1000 -m B -w > analysis/Ns25_pa_0.002_polyDFE.BIGRANGE.noDiv.out &

../../../poly_dfe/polyDFE/polyDFE -d analysis/Ns10_pa_0.005_polyDFE.txt -r ../../../poly_dfe/polyDFE/input/range_model_BandC.txt 1000 -m B > analysis/Ns10_pa_0.005_polyDFE.BIGRANGE.withDiv.out &
../../../poly_dfe/polyDFE/polyDFE -d analysis/Ns10_pa_0.005_polyDFE.txt -r ../../../poly_dfe/polyDFE/input/range_model_BandC.txt 1000 -m B -w > analysis/Ns10_pa_0.005_polyDFE.BIGRANGE.noDiv.out &

../../../poly_dfe/polyDFE/polyDFE -d analysis/Ns5_pa_0.01_polyDFE.txt -r ../../../poly_dfe/polyDFE/input/range_model_BandC.txt 1000 -m B > analysis/Ns5_pa_0.01_polyDFE.BIGRANGE.withDiv.out &
../../../poly_dfe/polyDFE/polyDFE -d analysis/Ns5_pa_0.01_polyDFE.txt -r ../../../poly_dfe/polyDFE/input/range_model_BandC.txt 1000 -m B -w > analysis/Ns5_pa_0.01_polyDFE.BIGRANGE.noDiv.out &
