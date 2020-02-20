# Here is where the initial conditions of the electron are defined
# This filename is the input parameter of the eTracks.py file

r_0 = .701142      # initial radius in c/w_p
pr_0 = 0.0            # initial momentum in m_e c
pz_0 = 0.0
xi_0 = 6.94922     # initial xi-position in c/w_p
SHModel = False    # True to use SHM model, False to use OSIRIS fields
track = 'min'
vr_0 = -1         # input arbitrary num when momenta are < 1 (takes care of gamma factor)
vz_0 = -1
