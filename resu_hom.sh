
# files 
vfil=(*.vtk)
nfil=${#vfil[*]}

for ((k=0 ; $nfil - $k ; k++))
do 

fil=${vfil[$k]}
./run_FFT $fil 22.393

done

