# Here is where the initial conditions of the electron are defined
# This filename is the input parameter of the eTracks.py file

r_0 = 0.779792     # 0.767132 #initial radius in c/w_p
vr_0 = (r_0 - 0.776168)/(837.322 - 837.317)
vz_0 = -1 # set arbitrarily when p < 1
pr_0 = 1.10027           # initial momentum in m_e c
pz_0 = -0.0445353
xi_0 = 6.47308    # initial xi-position in c/w_p
SHModel = False   # True to use SHM model, False to use OSIRIS fields
track = 'max'     # name of track if applicable
