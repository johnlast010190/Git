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
    (c) 2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "limitVelocity.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(limitVelocity, 0);
    addToRunTimeSelectionTable
    (
        option,
        limitVelocity,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::limitVelocity::limitVelocity
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    UName_(coeffs_.lookupOrDefault<word>("U", "U")),
    max_(coeffs_.lookup<scalar>("max"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::limitVelocity::sourceFields(wordList& fieldNames)
{
    fieldNames.setSize(1, UName_);
}


bool Foam::fv::limitVelocity::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        max_ = coeffs_.lookup<scalar>("max");

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::fv::limitVelocity::correct(volVectorField& U)
{
    const scalar maxSqrU = sqr(max_);

    vectorField& Uif = U.primitiveFieldRef();

    forAll(cells_, i)
    {
        const label celli = cells_[i];

        const scalar magSqrUi = magSqr(Uif[celli]);

        if (magSqrUi > maxSqrU)
        {
            Uif[celli] *= maxSqrU/magSqrUi;
        }
    }

    // handle boundaries in the case of 'all'
    if (selectionMode_ == smAll)
    {
        volVectorField::Boundary& Ubf = U.boundaryFieldRef();

        forAll(Ubf, patchi)
        {
            fvPatchVectorField& Up = Ubf[patchi];

            if (!Up.fixesValue())
            {
                forAll(Up, facei)
                {
                    const scalar magSqrUi = magSqr(Up[facei]);

                    if (magSqrUi > maxSqrU)
                    {
                        Up[facei] *= maxSqrU/magSqrUi;
                    }
                }
            }
        }
    }
}


// ************************************************************************* //
