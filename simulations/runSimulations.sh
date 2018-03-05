for i in $(ls)
do 
	mkdir $i/Rho0.009
	mkdir $i/Rho0.0045
	mkdir $i/Rho0.0009
	nice -n 10 parallel "slim $i/configs/Template.Rho0.009.txt | gzip > $i/Rho0.009/Rho0.009.{}.out" ::: $( seq 1 1000 )
	nice -n 10 parallel "slim $i/configs/Template.Rho0.0045.txt | gzip > $i/Rho0.0045/Rho0.0045.{}.out" ::: $( seq 1 1000 )
	nice -n 10 parallel "slim $i/configs/Template.Rho0.0009.txt | gzip > $i/Rho0.0009/Rho0.0009.{}.out" ::: $( seq 1 1000 )
done

