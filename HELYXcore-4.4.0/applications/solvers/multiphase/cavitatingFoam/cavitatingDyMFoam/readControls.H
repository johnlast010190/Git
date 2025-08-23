#include "cfdTools/general/include/readTimeControls.H"

maxAcousticCo = runTime.controlDict().lookup<scalar>("maxAcousticCo");
correctPhi = pimple.dict().lookupOrDefault<Switch>("correctPhi", true);
