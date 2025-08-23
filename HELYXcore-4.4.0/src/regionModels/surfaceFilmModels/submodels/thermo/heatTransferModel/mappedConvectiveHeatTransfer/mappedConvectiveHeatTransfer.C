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

\*---------------------------------------------------------------------------*/

#include "surfaceFilmModels/submodels/thermo/heatTransferModel/mappedConvectiveHeatTransfer/mappedConvectiveHeatTransfer.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "surfaceFilmModels/kinematicSingleLayer/kinematicSingleLayer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(mappedConvectiveHeatTransfer, 0);
addToRunTimeSelectionTable
(
    heatTransferModel,
    mappedConvectiveHeatTransfer,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mappedConvectiveHeatTransfer::mappedConvectiveHeatTransfer
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    heatTransferModel(film),
    htcConvPrimary_
    (
        IOobject
        (
            "htcConv",
            film.time().timeName(),
            film.primaryMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        film.primaryMesh()
    ),
    htcConvFilm_
    (
        IOobject
        (
            htcConvPrimary_.name(),
            film.time().timeName(),
            film.regionMesh()
        ),
        film.regionMesh(),
        dimensionedScalar(dimMass/pow3(dimTime)/dimTemperature, 0)
    )
{
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

mappedConvectiveHeatTransfer::~mappedConvectiveHeatTransfer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void mappedConvectiveHeatTransfer::correct()
{
    // Update the primary-side convective heat transfer coefficient
    htcConvPrimary_.correctBoundaryConditions();

    // Map the primary-side convective heat transfer coefficient
    // to the region internal field
    film().toRegion(htcConvFilm_, htcConvPrimary_.boundaryField());
}


tmp<volScalarField::Internal> mappedConvectiveHeatTransfer::h() const
{
    return htcConvFilm_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
