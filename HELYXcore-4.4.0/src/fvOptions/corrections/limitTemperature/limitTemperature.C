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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "limitTemperature.H"
#include "fvMesh/fvMesh.H"
#include "basicThermo/basicThermo.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(limitTemperature, 0);
    addToRunTimeSelectionTable
    (
        option,
        limitTemperature,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::limitTemperature::limitTemperature
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    TName_(coeffs_.lookupOrDefault<word>("field", "T")),
    Tmin_(coeffs_.lookup<scalar>("min")),
    Tmax_(coeffs_.lookup<scalar>("max"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::limitTemperature::sourceFields(wordList& fieldNames)
{
    // Set the field name to that of the energy field from which
    // the temperature is obtained
    if (obr_.foundObject<basicThermo>(basicThermo::dictName))
    {
        const basicThermo& thermo =
            obr_.lookupObject<basicThermo>(basicThermo::dictName);

        fieldNames.setSize(1, thermo.heT().name());
    }
    else
    {
        fieldNames.setSize(1, TName_);
    }
}


bool Foam::fv::limitTemperature::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        Tmin_ = coeffs_.lookup<scalar>("min");
        Tmax_ = coeffs_.lookup<scalar>("max");

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::fv::limitTemperature::correct(volScalarField& he)
{
    scalarField Tmin(cells_.size(), Tmin_);
    scalarField Tmax(cells_.size(), Tmax_);

    if (obr_.foundObject<basicThermo>(basicThermo::dictName))
    {
        const basicThermo& thermo =
            obr_.lookupObject<basicThermo>(basicThermo::dictName);
        const scalarField pSubSet(thermo.p(), cells_);
        scalarField heMin(thermo.he(pSubSet, Tmin, cells_));
        scalarField heMax(thermo.he(pSubSet, Tmax, cells_));

        scalarField& hec = he.primitiveFieldRef();

        forAll(cells_, i)
        {
            label celli = cells_[i];
            hec[celli] = max(min(hec[celli], heMax[i]), heMin[i]);
        }

        // handle boundaries in the case of 'all'
        if (selectionMode_ == smAll)
        {
            volScalarField::Boundary& bf = he.boundaryFieldRef();

            forAll(bf, patchi)
            {
                fvPatchScalarField& hep = bf[patchi];

                if (!hep.fixesValue())
                {
                    scalarField Tminp(hep.size(), Tmin_);
                    scalarField Tmaxp(hep.size(), Tmax_);

                    scalarField heMinp(thermo.he(Tminp, patchi));
                    scalarField heMaxp(thermo.he(Tmaxp, patchi));

                    forAll(hep, facei)
                    {
                        hep[facei] =
                            max(min(hep[facei], heMaxp[facei]), heMinp[facei]);
                    }
                }
            }
        }
    }
    else
    {
        scalarField& Tec = he.primitiveFieldRef();

        forAll(cells_, i)
        {
            label celli = cells_[i];
            Tec[celli] = max(min(Tec[celli], Tmax[i]), Tmin[i]);
        }

        // handle boundaries in the case of 'all'
        if (selectionMode_ == smAll)
        {
            volScalarField::Boundary& bf = he.boundaryFieldRef();

            forAll(bf, patchi)
            {
                fvPatchScalarField& Tep = bf[patchi];

                if (!Tep.fixesValue())
                {
                    scalarField Tminp(he.size(), Tmin_);
                    scalarField Tmaxp(he.size(), Tmax_);

                    forAll(Tep, facei)
                    {
                        Tep[facei] =
                            max(min(Tep[facei], Tmaxp[facei]), Tminp[facei]);
                    }
                }
            }
        }
    }
}


// ************************************************************************* //
