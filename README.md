# PWFA-eTracks

This code was developed to simulate electrons injected into the wake of an electron-driven plasma wakefield accelerator. It does so by taking previously simulated fields from OSIRIS and propagating an electron through these fields. This allows us to quickly estimate the trajectories of individual electrons without running OSIRIS repeatedly, which would be computationally expensive.

The master branch contains the most up-to-date code for tracking single electron trajectories

### Running the Simulation
To run without a pre-defined initial condition file use
```
python3 eTracks.py
```
and you will be prompted to input the initial conditions in normalized units.

You can also add an initialization file in the `input` folder, such as the file `example.py`, which would be input to the simulation using 
```
python3 eTracks.py inputs.example
```
### Requirements
This simulation requires the python packages h5py, importlib, numpy, and matplotlib.

### Contact
Contact Audrey Farrell (audreyclairefarrell@gmail.com) for questions about this code. Source code can be found at https://github.com/acfarrell/PWFA-eTracks
