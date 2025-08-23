/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
|    o     o     |  ENGYS Ltd. <http://engys.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of HELYXcore.
    HELYXcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    HELYXcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HELYXcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with HELYXcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2021-2025 Engys Ltd.

Class
    Foam::materialModels

Description
    Interface classes.

SourceFiles
    materialModels.C

\*---------------------------------------------------------------------------*/


#include "fields/volFields/volFields.H"
#include "baseModels.H"
#include "materialModels/materialMacros.H"

#ifndef materialModels_H
#define materialModels_H
    defineScalarMaterialModel(W);
    defineScalarMaterialModel(Y);
    defineScalarMaterialModel(R);
    defineScalarMaterialModel(limit);
    defineScalarMaterialModel(rho);
    defineScalarMaterialModel(hContribution);
    defineScalarMaterialModel(CpContribution);
    defineScalarMaterialModel(eContribution);
    defineScalarMaterialModel(CvContribution);
    defineScalarMaterialModel(spContribution);
    defineScalarMaterialModel(svContribution);
    defineScalarMaterialModel(psi);
    defineScalarMaterialModel(Z);
    defineScalarMaterialModel(CpMCv);
    defineScalarMaterialModel(hf);
    defineScalarMaterialModel(a);
    defineScalarMaterialModel(Cp);
    defineScalarMaterialModel(Cpv);
    defineScalarMaterialModel(Cv);
    defineScalarMaterialModel(ea);
    defineScalarMaterialModel(es);
    defineScalarMaterialModel(gamma);
    defineScalarMaterialModel(g);
    defineScalarMaterialModel(ha);
    defineScalarMaterialModel(he);
    defineScalarMaterialModel(hs);
    defineScalarMaterialModel(s);
    defineScalarMaterialModel(strainRate);
    defineScalarMaterialModel(mu);
    defineScalarMaterialModel(kappa);
    defineVectorMaterialModel(vKappa);
    defineScalarMaterialModel(weightedAverageMixture);
    defineVectorMaterialModel(vectorWeightedAverageMixture);
    defineTensorMaterialModel(tensorWeightedAverageMixture);
    defineScalarMaterialModel(clippedWeightedAverageMixture);
    defineVectorMaterialModel(vectorClippedWeightedAverageMixture);
    defineTensorMaterialModel(tensorClippedWeightedAverageMixture);
    defineScalarMaterialModel(harmonicMixture);
    defineScalarMaterialModel(egrReactantMixture);
    defineScalarMaterialModel(egrProductMixture);
    defineScalarMaterialModel(egrThermoMixture);
    defineScalarMaterialModel(homogeneousReactantMixture);
    defineScalarMaterialModel(homogeneousProductMixture);
    defineScalarMaterialModel(homogeneousThermoMixture);
    defineScalarMaterialModel(inhomogeneousReactantMixture);
    defineScalarMaterialModel(inhomogeneousProductMixture);
    defineScalarMaterialModel(inhomogeneousThermoMixture);
    defineScalarMaterialModel(veryInhomogeneousReactantMixture);
    defineScalarMaterialModel(veryInhomogeneousProductMixture);
    defineScalarMaterialModel(veryInhomogeneousThermoMixture);
    defineScalarMaterialModel(T);
    defineScalarMaterialModel(Tha);
    defineScalarMaterialModel(Ths);
    defineScalarMaterialModel(Tea);
    defineScalarMaterialModel(Tes);
    defineScalarMaterialModel(p);
    defineScalarMaterialModel(dCpdT);
    defineScalarMaterialModel(alphav);

    // Liquid properties
    defineScalarMaterialModel(pv);
    defineScalarMaterialModel(hl);
    defineScalarMaterialModel(Cpg);
    defineScalarMaterialModel(B);
    defineScalarMaterialModel(kappag);
    defineScalarMaterialModel(mug);
    defineScalarMaterialModel(sigma);
    defineScalarMaterialModel(D);
    defineScalarMaterialModel(pvInvert);
    defineScalarMaterialModel(c0);

    // Sturation models
    defineScalarMaterialModel(pSat);
    defineScalarMaterialModel(pSatPrime);
    defineScalarMaterialModel(lnPSat);
    defineScalarMaterialModel(Tsat);

    // Solid properties
    defineScalarMaterialModel(E);
    defineScalarMaterialModel(nu);

#endif

// ************************************************************************* //
