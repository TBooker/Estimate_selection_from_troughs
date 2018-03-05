for i in $(seq 1 1000)
	do
	echo $i
#	mv Rho0.0045.$i.out Rho0.0045.$i.out.gz
#	mv Exp.Nes100_dDFE.Exponential_div10.Rho0.0009.$i.out.gz Rho0.0009.$i.out.gz
	mv Exp.Nes10_dDFE.Exponential.Rho0.0045.$i.out.gz Rho0.0045.$i.out.gz
	done
