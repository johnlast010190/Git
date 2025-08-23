nCorr = stressControl.lookupOrDefault<int>("nCorrectors", 1);
convergenceTolerance = stressControl.lookup<scalar>("D");
compactNormalStress = stressControl.lookup<Switch>("compactNormalStress");
