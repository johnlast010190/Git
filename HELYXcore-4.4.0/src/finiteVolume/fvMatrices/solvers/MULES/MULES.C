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
    (c) 2016 OpenCFD Ltd.
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fvMatrices/solvers/MULES/MULES.H"
#if !defined( WIN32 ) && !defined( WIN64 )
#include "global/profiling/profiling.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::MULES::explicitSolve
(
    volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiPsi,
    const scalar psiMax,
    const scalar psiMin
)
{
    const dictionary& MULEScontrols =
        psi.mesh().solution().solverDict(psi.name());
    explicitSolve(psi, phi, phiPsi, psiMax, psiMin, MULEScontrols);
}


void Foam::MULES::explicitSolve
(
    volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiPsi,
    const scalar psiMax,
    const scalar psiMin,
    const dictionary& MULEScontrols
)
{
    #if !defined( WIN32 ) && !defined( WIN64 )
    addProfiling(solve, "MULES::explicitSolve");
    #endif

    explicitSolve
    (
        geometricOneField(),
        psi,
        phi,
        phiPsi,
        zeroField(), zeroField(),
        psiMax, psiMin,
        MULEScontrols
    );
}


void Foam::MULES::limitSum(UPtrList<scalarField>& phiPsiCorrs)
{
    forAll(phiPsiCorrs[0], facei)
    {
        scalar sumPos = 0;
        scalar sumNeg = 0;

        for (int phasei=0; phasei<phiPsiCorrs.size(); phasei++)
        {
            if (phiPsiCorrs[phasei][facei] > 0)
            {
                sumPos += phiPsiCorrs[phasei][facei];
            }
            else
            {
                sumNeg += phiPsiCorrs[phasei][facei];
            }
        }

        scalar sum = sumPos + sumNeg;

        if (sum > 0 && sumPos > VSMALL)
        {
            scalar lambda = -sumNeg/sumPos;

            for (int phasei=0; phasei<phiPsiCorrs.size(); phasei++)
            {
                if (phiPsiCorrs[phasei][facei] > 0)
                {
                    phiPsiCorrs[phasei][facei] *= lambda;
                }
            }
        }
        else if (sum < 0 && sumNeg < -VSMALL)
        {
            scalar lambda = -sumPos/sumNeg;

            for (int phasei=0; phasei<phiPsiCorrs.size(); phasei++)
            {
                if (phiPsiCorrs[phasei][facei] < 0)
                {
                    phiPsiCorrs[phasei][facei] *= lambda;
                }
            }
        }
    }
}


// ************************************************************************* //
