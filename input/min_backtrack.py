# Here is where the initial conditions of the electron are defined
# This filename is the input parameter of the eTracks.py file

r_0 = 0.201935      # initial radius in c/w_p
pr_0 = 0.000332963            # initial momentum in m_e c
pz_0 = 97.0326
xi_0 = 1.32031     # initial xi-position in c/w_p
SHModel = False    # True to use SHM model, False to use OSIRIS fields
track = 'min'
vr_0 = -1         # input arbitrary num when momenta are < 1 (takes care of gamma factor)
vz_0 = -1
