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
    (c) 2018 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"
#include "accelerationSource.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(accelerationSource, 0);
    addToRunTimeSelectionTable(option, accelerationSource, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::accelerationSource::accelerationSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    velOrAccel_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::accelerationSource::sourceFields(wordList& fieldNames)
{
    fieldNames = wordList(1, coeffs_.lookupOrDefault<word>("U", "U"));
}


void Foam::fv::accelerationSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    add(geometricOneField(), eqn, fieldi);
}


void Foam::fv::accelerationSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    add(rho, eqn, fieldi);
}


void Foam::fv::accelerationSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    add((alpha*rho)(), eqn, fieldi);
}


void Foam::fv::accelerationSource::addSup
(
    fvBlockMatrix<vector>& eqn,
    const label fieldi
)
{
    add(geometricOneField(), eqn, fieldi);
}


void Foam::fv::accelerationSource::addSup
(
    const volScalarField& rho,
    fvBlockMatrix<vector>& eqn,
    const label fieldi
)
{
    add(rho, eqn, fieldi);
}


bool Foam::fv::accelerationSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        if (dict.found("velocity"))
        {
            velOrAccel_ = Function1<vector>::New("velocity", dict);
        }
        else if (dict.found("acceleration"))
        {
            velOrAccel_ = Function1<vector>::New("acceleration", dict);
        }
        else
        {
            FatalErrorInFunction
                << "velocity or acceleration entry is not defined."
                << abort(FatalError);
        }

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
