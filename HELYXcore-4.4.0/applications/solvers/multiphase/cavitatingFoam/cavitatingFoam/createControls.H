#include "cfdTools/general/include/createTimeControls.H"

scalar maxAcousticCo
(
    runTime.controlDict().lookup<scalar>("maxAcousticCo")
);
