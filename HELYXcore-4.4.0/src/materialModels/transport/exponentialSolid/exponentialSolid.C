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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2022-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "exponentialSolid.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(exponentialSolid, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        exponentialSolid,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::exponentialSolid::exponentialSolid
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
:
    materialModel(obr, dict, phaseName, specieName, name)
{
    exponentialSolid::read();
}


Foam::autoPtr<Foam::exponentialSolid>
Foam::exponentialSolid::clone() const
{
    return autoPtr<exponentialSolid>(new exponentialSolid(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::baseModels<Foam::scalar>* Foam::exponentialSolid::castScalarModel
(
    const word& modelName
)
{
    castMaterial(modelName, kappa)

    return nullptr;
}


Foam::baseModels<Foam::vector>* Foam::exponentialSolid::castVectorModel
(
    const word& modelName
)
{
    castMaterial(modelName, vKappa)

    return nullptr;
}


bool Foam::exponentialSolid::read()
{
    initialise("T");
    kappa0_ = dict_->lookup<scalar>("kappa0");
    n0_ = dict_->lookup<scalar>("n0");
    Tref_ = dict_->lookup<scalar>("Tref")+T_->offset().value();

    return true;
}


// ************************************************************************* //
