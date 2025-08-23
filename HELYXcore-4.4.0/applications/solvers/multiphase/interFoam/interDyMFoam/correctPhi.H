CorrectPhi
(
    U,
    phi,
    p_rgh,
    surfaceScalarField("rAUf", fvc::interpolate(rAU)),
    geometricZeroField(),
    pimple,
    pRefCell
);

#include "cfdTools/incompressible/continuityErrs.H"
