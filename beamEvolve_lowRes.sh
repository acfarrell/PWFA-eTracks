#!/bin/bash

#Reset initial input file to the first timestep (100)
cp input/quickPIC_100_lowRes.py input/quickPIC_lowRes.py
old_t=100
for t in {100..200}
do
  echo -ne "Simulating t = $t\r"
  # Edit initial file to update filenames and time step
  fname=slicexz_00000${t}.h5
  old_fname=slicexz_00000${old_t}.h5
  sed -i 's/t0 = '"$old_t"'/t0 = '"$t"'/g' input/quickPIC_lowRes.py
  sed -i 's/'"$old_fname"'/'"$fname"'/g' input/quickPIC_lowRes.py

  output_fname=quickPIC_$t
  old_output_fname=quickPIC_$old_t
  sed -i 's/'"$old_output_fname"'/'"$output_fname"'/g' input/quickPIC_lowRes.py
  old_t=$t
  python3 beamSim.py input.quickPIC_lowRes
done
