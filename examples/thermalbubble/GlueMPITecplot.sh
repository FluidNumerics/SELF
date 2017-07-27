#!/bin/bash

# Obtain the number of MPI processes
nP=4
#$( grep nProc runtime.params | awk '{gsub(/,$/,""); print $3}' )
iter0=$( grep iterInit runtime.params | awk '{gsub(/,$/,""); print $3}' )
dFreq=$( grep dumpFreq runtime.params | awk '{gsub(/,$/,""); print $3}' )
nT=$( grep nTimeSteps runtime.params | awk '{gsub(/,$/,""); print $3}' )
iterEnd=$( echo "${nT}+${iter0}" | bc )
nPu=$( echo "${nP}-1" | bc)

echo 'runtime.params : nProcessess =' ${nP}
echo 'runtime.params  : iterInit    = '${iter0}
echo 'runtimes.params : dumpFreq    =' ${dFreq}
echo 'runtimes.params : nTimeSteps  =' ${nT}
echo 'Last Iterate    : iterEnd     =' ${iterEnd}

# Loop over model iterates
for i in $( seq -f "%010g" $iter0 $dFreq $iterEnd )
do 
   tname='mastertemp.tec'
   fpname='Euler.0000.'$i'.tec'
   cat $fpname > $tname
   #rm $fpname
   for j in $( seq -f "%04g" 1 $nPu )
   do
      fpname='Euler.'$j'.'$i'.tec'
      echo $fpname
      # Temporarily copy the "master" file
      cat $tname > temp1.tec
      # Copy all except the first line of the tec file
      cat $fpname | sed 1d > temp2.tec
      # Concatenate the master file with the process file
      cat temp1.tec temp2.tec > $tname
      #rm $fpname
   done
   fname='Euler.'$i'.tec'
   cat mastertemp.tec > $fname
done

rm temp1.tec temp2.tec mastertemp.tec
