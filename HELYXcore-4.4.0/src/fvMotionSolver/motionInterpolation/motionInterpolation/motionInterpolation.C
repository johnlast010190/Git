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
    (c) 2015 OpenCFD Ltd.
    (c) 2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "motionInterpolation/motionInterpolation/motionInterpolation.H"
#include "interpolation/volPointInterpolation/volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(motionInterpolation, 0);

    defineRunTimeSelectionTable(motionInterpolation, Istream);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::motionInterpolation::motionInterpolation
(
    const fvMesh& mesh
)
:
    mesh_(mesh)
{}


Foam::motionInterpolation::motionInterpolation
(
    const fvMesh& mesh,
    Istream& entry
)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::motionInterpolation>
Foam::motionInterpolation::New(const fvMesh& mesh)
{
    return autoPtr<motionInterpolation>(new motionInterpolation(mesh));
}


Foam::autoPtr<Foam::motionInterpolation>
Foam::motionInterpolation::New(const fvMesh& mesh, Istream& entry)
{
    const word type(entry);

    Info<< "Selecting motion interpolation: " << type << endl;

    const auto ctor =
        ctorTableLookup
        (
            "interpolation type",
            IstreamConstructorTable_(),
            type
        );
    return autoPtr<motionInterpolation>(ctor(mesh, entry));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::motionInterpolation::~motionInterpolation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::motionInterpolation::interpolate
(
    const volScalarField& cellDisplacement,
    pointScalarField& pointDisplacement
) const
{
    volPointInterpolation::New(mesh()).interpolate
    (
        cellDisplacement,
        pointDisplacement
    );
}


void Foam::motionInterpolation::interpolate
(
    const volVectorField& cellDisplacement,
    pointVectorField& pointDisplacement
) const
{
    volPointInterpolation::New(mesh()).interpolate
    (
        cellDisplacement,
        pointDisplacement
    );
}


// ************************************************************************* //
