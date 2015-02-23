
# This script is a rewrite of the matlab twins code. 
# The purpose of the rewrite is to make a self-contained script
# that doesn't require matlab to be run and to make the calculation
# simpler. 

### PHYSICAL PARAMETERS ###

lambda0 = 8E-9 #[m]  the lambdas are the mean free paths 
lambda1 = 30E-9 #[m]
ncarrier_0 = 2.8E16 #[cm^-3] 'undoped' carrier density
(mu_min, mu_max) = (68.5, 1414) # mobility limits for Si, more info
# at http://www.cleanroom.byu.edu/ResistivityCal.phtml 

