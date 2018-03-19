python ~/project/15.Estimate_Selection_From_Trough/simulations/run_polyDFE/bin/generateInputFiles.py -i $1 -b $2 -o polyDFE.config.$2.txt -s m6 --polyDFE
/home/booker/project/15.Estimate_Selection_From_Trough/poly_dfe/polyDFE/polyDFE -d polyDFE.config.$2.txt -m C > polyDFE.output.$2.txt
