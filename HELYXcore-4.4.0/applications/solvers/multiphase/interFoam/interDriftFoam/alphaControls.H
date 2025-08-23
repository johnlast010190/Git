const dictionary& alphaControls = mesh.solution().solverDict(alpha1.name());

label nAlphaCorr(alphaControls.lookup<label>("nAlphaCorr"));

label nAlphaSubCycles(alphaControls.lookup<label>("nAlphaSubCycles"));

scalar blend(alphaControls.lookupOrDefault<scalar>("blend", scalar(0.99)));
