for ((i=0;i<=100;i++));do file=`printf "micro_%04d" $i`;./write_vtk.exe $file;done
