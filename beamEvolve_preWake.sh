#!/bin/bash

#Reset initial input file to the first timestep (100)
cp input/quickPIC_100.py input/quickPIC.py
old_t=100
for t in {0..99}
do
  echo -ne "Simulating t = $t\r"
  # Edit initial file to update filenames and time step
  fname=$(printf "slicexz_00000%03d" ${t}).h5
  old_fname=$(printf "slicexz_00000%03d" ${old_t}).h5
  echo -ne "Reading File $fname\r"
  sed -i 's/t0 = '"$old_t"'/t0 = '"$t"'/g' input/quickPIC.py
  sed -i 's/'"$old_fname"'/'"$fname"'/g' input/quickPIC.py

  output_fname=quickPIC_$t
  old_output_fname=quickPIC_$old_t
  sed -i 's/'"$old_output_fname"'/'"$output_fname"'/g' input/quickPIC.py
  old_t=$t
  python3 beamSim.py input.quickPIC
done
