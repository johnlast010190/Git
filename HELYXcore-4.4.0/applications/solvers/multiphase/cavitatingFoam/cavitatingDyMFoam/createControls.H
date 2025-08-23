#include "cfdTools/general/include/createTimeControls.H"

scalar maxAcousticCo
(
    runTime.controlDict().lookup<scalar>("maxAcousticCo")
);

bool correctPhi
(
    pimple.dict().lookupOrDefault<Switch>("correctPhi", true)
);
