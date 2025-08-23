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
    (c) 2016 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "meanSurfaceIntensity.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "basicThermo/basicThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fieldInitializations
{
    defineTypeNameAndDebug(meanSurfaceIntensity, 0);
    addToRunTimeSelectionTable(fieldInit, meanSurfaceIntensity, initMethod);
}
}

// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldInitializations::meanSurfaceIntensity::meanSurfaceIntensity
(
    const fvMesh& mesh,
    const fvSolutionRegistry& localDb,
    const dictionary& fieldDict,
    const word& fN
)
:
    fieldInit(mesh, localDb, fieldDict, fN)
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fieldInitializations::meanSurfaceIntensity::correct()
{
    //if already initialised stop
    if (initialised())
    {
        return;
    }

    fieldInit::initMsg();

    word TName = initDict().lookupOrDefault<word>("T", "T");
    scalar sigmaSB = initDict().lookupOrDefault<scalar>("sigmaSB", 5.6704e-8);
    scalar emissivity = initDict().lookupOrDefault<scalar>("emissivity", 1);

    volScalarField& T =
        const_cast<volScalarField&>
        (localDb().lookupObject<volScalarField>(TName));

    scalar TwMean = 0;
    scalar Awall = 0;

    forAll(T.boundaryField(), bI)
    {
        if (isA<wallFvPatch>(mesh().boundary()[bI]))
        {
            TwMean
                += gSum(T.boundaryField()[bI]
                * mesh().boundary()[bI].magSf());
            Awall += gSum(mesh().boundary()[bI].magSf());
        }
    }
    TwMean /= Awall;

    volScalarField& f
        = const_cast<volScalarField&>
        (localDb().lookupObject<volScalarField>(name()));

    f.primitiveFieldRef() =
        emissivity * sigmaSB
      * Foam::pow(TwMean+basicThermo::TRefIfFound(localDb()).value(), 4.0)
      / constant::mathematical::pi;

    // the field has been initialised
    initialised() = true;

    initCoupledBoundaries();
}



// ************************************************************************* //
