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
    (c) 2013-2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "submodels/MPPIC/TimeScaleModels/TimeScaleModel/TimeScaleModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(TimeScaleModel, 0);
    defineRunTimeSelectionTable(TimeScaleModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TimeScaleModel::TimeScaleModel
(
    const dictionary& dict
)
:
    alphaPacked_(dict.lookup<scalar>("alphaPacked")),
    e_(dict.lookup<scalar>("e"))
{
}


Foam::TimeScaleModel::TimeScaleModel
(
    const TimeScaleModel& cm
)
:
    alphaPacked_(cm.alphaPacked_),
    e_(cm.e_)
{
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::TimeScaleModel> Foam::TimeScaleModel::New
(
    const dictionary& dict
)
{
    word modelType(dict.lookup("type"));

    Info<< "Selecting time scale model " << modelType << endl;

    const auto ctor =
        ctorTableLookup
        (
            "time scale model type",
            dictionaryConstructorTable_(),
            modelType
        );
    return autoPtr<TimeScaleModel>(ctor(dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TimeScaleModel::~TimeScaleModel()
{}


// ************************************************************************* //
