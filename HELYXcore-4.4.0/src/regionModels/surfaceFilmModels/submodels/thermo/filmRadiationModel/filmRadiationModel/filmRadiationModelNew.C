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
    (c) 2011-2020 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "surfaceFilmModels/submodels/thermo/filmRadiationModel/filmRadiationModel/filmRadiationModel.H"
#include "surfaceFilmModels/submodels/thermo/filmRadiationModel/noRadiation/noRadiation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<radiationModel> radiationModel::New
(
    surfaceFilmRegionModel& model,
    const dictionary& dict
)
{
    if
    (
        !dict.found("radiationModel")
     && !dict.found(radiationModel::typeName)
    )
    {
        return autoPtr<radiationModel>(new noRadiation(model, dict));
    }

    const dictionary& radiationDict
    (
        dict.found("radiationModel")
      ? dict
      : dict.subDict(radiationModel::typeName)
    );

    const word modelType
    (
        dict.found("radiationModel")
      ? radiationDict.lookup("radiationModel")
      : radiationDict.lookup("model")
    );

    Info<< "    Selecting radiationModel " << modelType << endl;

    const auto ctor =
        ctorTableLookup
        (
            "radiationModel type",
            dictionaryConstructorTable_(),
            modelType
        );
    return autoPtr<radiationModel>(ctor(model, radiationDict));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
