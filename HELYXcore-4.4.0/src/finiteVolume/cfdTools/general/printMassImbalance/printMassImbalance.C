
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
    (c) 2018 Engys Ltd.

Class
    Foam::printMassImbalance

Description

SourceFiles
    printMassImbalance.H

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/printMassImbalance/printMassImbalance.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::printMassImbalance
(
    const word& regionName,
    const surfaceScalarField& phi
)
{

    scalar massIn = 0.0;
    scalar massOut = 0.0;
    //- take GIBs into account
    tmp<surfaceScalarField> tphi(fvc::applyFaceMask(phi));
    const surfaceScalarField::Boundary& bPhi = tphi().boundaryField();

    //BC in-out
    forAll(bPhi, pI)
    {
        if (!bPhi[pI].coupled())
        {
            const fvsPatchScalarField& phip = bPhi[pI];
            forAll(bPhi[pI], fI)
            {
                if (phip[fI] < 0.0) massIn -= phip[fI];
                else if (phip[fI] > 0.0) massOut += phip[fI];
            }
        }
    }
    reduce(
        std::tie(massIn, massOut),
        UniformParallelOp<sumOp<scalar>, 2>{}
    );
    if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        Info<< "Mass-imbalance [kg/s] (";
    }
    else if (phi.dimensions() == dimVelocity*dimArea)
    {
        Info<< "Mass-imbalance [m^3/s] (";
    }
    else
    {
        FatalErrorInFunction
            << "Wrong phi dimensions " << phi.dimensions()
            << exit(FatalError);
    }
    Info<< regionName << ")"
         << ": massIn = " << massIn
         << ", massOut = " << massOut
         << ", Delta = " << massIn-massOut
         << endl;
}

// ************************************************************************* //
