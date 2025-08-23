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
    (c) 2012-2019 OpenFOAM Foundation
    (c) 2010-2016 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "IDDESDelta2.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(IDDESDelta2, 0);
    addToRunTimeSelectionTable(LESdelta, IDDESDelta2, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::IDDESDelta2::calcDelta()
{
    const volScalarField& hmax = hmax_;
    const fvMesh& mesh = turbulenceModel_.mesh();

    // Wall-normal vectors
    const volVectorField& n = wallDist::New(mesh).n();

    tmp<volScalarField> tHwallNormal
    (
        volScalarField::New
        (
            "HwallNormal",
            mesh,
            dimensionedScalar(dimLength, 0)
        )
    );

    scalarField& HwallNormal = tHwallNormal.ref().primitiveFieldRef();

    const cellList& cells = mesh.cells();
    const vectorField& Sf(mesh.faceAreas());
    const scalarField& Vol(mesh.V());

    forAll(cells, cellI)
    {
        scalar Snorm = 0;
        const labelList& cFaces = cells[cellI];
        const vector nCell = n[cellI];
        forAll(cFaces, cFaceI)
        {
            label faceI = cFaces[cFaceI];
            Snorm += mag(Sf[faceI] & nCell);
        }
        HwallNormal[cellI] = 2*Vol[cellI] / Snorm;
    }


    label nD = mesh.nGeometricD();

    if (nD == 2)
    {
        WarningInFunction
            << "Case is 2D, LES is not strictly applicable" << nl
            << endl;
    }
    else if (nD != 3)
    {
        FatalErrorInFunction
            << "Case must be either 2D or 3D" << exit(FatalError);
    }

    delta_.primitiveFieldRef() =
        min
        (
            max
            (
                max
                (
                    Cw_*wallDist::New(mesh).y(),
                    Cw_*hmax
                ),
                tHwallNormal()
            ),
            hmax
        );

    delta_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESModels::IDDESDelta2::IDDESDelta2
(
    const word& name,
    const turbulenceModel& turbulence,
    const dictionary& dict
)
:
    LESdelta(name, turbulence),
    hmax_
    (
        IOobject::groupName("hmax", turbulence.U().group()),
        turbulence,
        dict
    ),
    Cw_
    (
        dict.subDict(type() + "Coeffs").lookupOrDefault<scalar>("Cw", 0.15)
    )
{
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::IDDESDelta2::read(const dictionary& dict)
{
    const dictionary& coeffsDict(dict.subDict(type() + "Coeffs"));

    coeffsDict.readIfPresent<scalar>("Cw", Cw_);

    calcDelta();
}


void Foam::LESModels::IDDESDelta2::correct()
{
    if (turbulenceModel_.mesh().changing())
    {
        calcDelta();
        hmax_.correct();
    }
}


// ************************************************************************* //
