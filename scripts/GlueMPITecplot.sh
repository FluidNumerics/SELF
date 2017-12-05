#!/bin/bash

# Obtain the number of MPI processes
nP=8
nPu=$( echo "${nP}-1" | bc)

echo 'runtime.params : nProcessess =' ${nP}

# Loop over model iterates
for file in $( ls State.0000.*.tec )
do 
   tname='mastertemp.tec'
   cat $file > $tname

   for j in $( seq -f "%04g" 1 $nPu )
   do
      iter=$( echo $file | awk -F '.' '{print $3}' )
      fpname='State.'$j'.'$iter'.tec'
      echo $fpname
      # Temporarily copy the "master" file
      cat $tname > temp1.tec
      # Copy all except the first line of the tec file
      cat $fpname | sed 1d > temp2.tec
      # Concatenate the master file with the process file
      cat temp1.tec temp2.tec > $tname
   #   rm $fpname
   done
   fname='State.'$iter'.tec'
   cat mastertemp.tec > $fname
done

rm temp1.tec temp2.tec mastertemp.tec
