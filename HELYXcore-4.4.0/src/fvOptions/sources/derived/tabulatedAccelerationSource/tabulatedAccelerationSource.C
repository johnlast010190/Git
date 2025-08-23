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
    (c) 2015-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "tabulatedAccelerationSource.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(tabulatedAccelerationSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        tabulatedAccelerationSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::tabulatedAccelerationSource::tabulatedAccelerationSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    option(name, modelType, dict, obr),
    motion_(coeffs_, mesh_.time()),
    UName_(coeffs_.lookupOrDefault<word>("U", "U")),
    g0_("g0", dimAcceleration, Zero)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::tabulatedAccelerationSource::initialise()
{
    if (mesh().foundObject<uniformDimensionedVectorField>("g"))
    {
        g0_ = obr_.lookupObject<uniformDimensionedVectorField>("g");
    }

    return true;
}


void Foam::fv::tabulatedAccelerationSource::sourceFields
(
    wordList& fieldNames
)
{
    fieldNames.setSize(1, UName_);
}


void Foam::fv::tabulatedAccelerationSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    addSup<geometricOneField>(geometricOneField(), eqn, fieldi);
}


void Foam::fv::tabulatedAccelerationSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    addSup<volScalarField>(rho, eqn, fieldi);
}


bool Foam::fv::tabulatedAccelerationSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        return motion_.read(coeffs_);
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
