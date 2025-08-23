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
    (c) 2015-2019 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "absorptionEmissionModels/multiBandZoneAbsorptionEmission/multiBandZoneAbsorptionEmission.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{
    defineTypeNameAndDebug(multiBandZoneAbsorptionEmission,
    0);

    addToRunTimeSelectionTable
    (
        absorptionEmissionModel,
        multiBandZoneAbsorptionEmission,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::multiBandZoneAbsorptionEmission::
multiBandZoneAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& typeNameDerived
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeNameDerived + "Coeffs")),
    absCoeffs_(maxBands_),
    emiCoeffs_(maxBands_),
    nBands_(0),
    zoneAbsorptivity_(),
    zoneEmissivity_(),
    zoneIds_()
{
    coeffsDict_.readEntry("absorptivity", absCoeffs_);
    coeffsDict_.readEntry("emissivity", emiCoeffs_);
    nBands_ = absCoeffs_.size();

    const dictionary& zoneDict = coeffsDict_.subDict("zones");

    zoneDict.readEntry("absorptivity", zoneAbsorptivity_);
    zoneDict.readEntry("emissivity", zoneEmissivity_);

    zoneIds_.resize(zoneAbsorptivity_.size(), -1);

    label numZones = 0;
    forAllConstIters(zoneAbsorptivity_, iter)
    {
        label zoneID = mesh.cellZones().findZoneID(iter.key());
        if (zoneID == -1)
        {
            FatalErrorInFunction
                << "Cannot find cellZone " << iter.key() << endl
                << "Valid cellZones are " << mesh.cellZones().names()
                << exit(FatalError);
        }
        zoneIds_[numZones] = zoneID;
        ++numZones;
    }
    // zoneIds_.resize(numZones);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::multiBandZoneAbsorptionEmission::
~multiBandZoneAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::multiBandZoneAbsorptionEmission::
aCont
(
    const label bandI
) const
{
    auto ta = volScalarField::New
    (
        "aCont",
        mesh(),
        dimensionedScalar("a", dimless/dimLength, absCoeffs_[bandI])
    );
    scalarField& a = ta.ref().primitiveFieldRef();

    for (const label zonei : zoneIds_)
    {
        const cellZone& zn = mesh().cellZones()[zonei];
        const auto iter = zoneAbsorptivity_.cfind(zn.name());

        if (iter.found())  // Check is redundant (cannot fail)
        {
            const scalarList& absorb = iter.object();

            UIndirectList<scalar>(a, zn) = absorb[bandI];
        }
    }

    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::multiBandZoneAbsorptionEmission::
eCont
(
    const label bandI
) const
{
    auto te = volScalarField::New
    (
        "eCont",
        mesh(),
        dimensionedScalar("e", dimless/dimLength, emiCoeffs_[bandI])
    );
    scalarField& e = te.ref().primitiveFieldRef();

    for (const label zonei : zoneIds_)
    {
        const cellZone& zn = mesh().cellZones()[zonei];
        const auto iter = zoneEmissivity_.cfind(zn.name());

        if (iter.found())  // Check is redundant (cannot fail)
        {
            const scalarList& emit = iter.object();

            UIndirectList<scalar>(e, zn) = emit[bandI];
        }
    }

    return te;
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::multiBandZoneAbsorptionEmission::
ECont
(
    const label bandI
) const
{
    return volScalarField::New
    (
        "ECont",
        mesh(),
        dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), Zero)
    );
}


// ************************************************************************* //
