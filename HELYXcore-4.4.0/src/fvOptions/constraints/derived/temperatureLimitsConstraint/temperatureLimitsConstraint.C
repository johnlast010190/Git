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
    (c) 2012-2013 OpenFOAM Foundation
    (c) 2024 Engys Ltd

\*----------------------------------------------------------------------------*/

#include "temperatureLimitsConstraint.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "basicThermo/basicThermo.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(temperatureLimitsConstraint, 0);
    addToRunTimeSelectionTable
    (
        option,
        temperatureLimitsConstraint,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::temperatureLimitsConstraint::temperatureLimitsConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    Tmin_(coeffs_.lookup<scalar>("Tmin")),
    Tmax_(coeffs_.lookup<scalar>("Tmax"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::temperatureLimitsConstraint::sourceFields(wordList& fieldNames)
{
    const basicThermo& thermo =
        obr_.lookupObject<basicThermo>(basicThermo::dictName);
    fieldNames.setSize(1, thermo.heT().name());
}


bool Foam::fv::temperatureLimitsConstraint::alwaysApply() const
{
    return true;
}


void Foam::fv::temperatureLimitsConstraint::correct(volScalarField& heT)
{
    const basicThermo& thermo =
        obr_.lookupObject<basicThermo>(basicThermo::dictName);

    if (heT.name() == thermo.he().name())
    {
        scalarField Tmin(cells_.size(), Tmin_);
        scalarField Tmax(cells_.size(), Tmax_);
        const scalarField pSubSet(thermo.p(), cells_);
        scalarField heMin(thermo.he(pSubSet, Tmin, cells_));
        scalarField heMax(thermo.he(pSubSet, Tmax, cells_));

        scalarField& hec = heT.ref();

        forAll(cells_, i)
        {
            label cellI = cells_[i];
            hec[cellI] = max(min(hec[cellI], heMax[i]), heMin[i]);
        }

        // handle boundaries in the case of 'all'
        if (selectionMode_ == smAll)
        {
            volScalarField::Boundary& bf = heT.boundaryFieldRef();

            forAll(bf, patchI)
            {
                fvPatchScalarField& hep = bf[patchI];
                if (hep.fixesValue())
                {
                    // not over-riding fixed conditions
                    continue;
                }

                scalarField Tminp(hep.size(), Tmin_);
                scalarField Tmaxp(hep.size(), Tmax_);

                scalarField heMinp(thermo.he(Tminp, patchI));
                scalarField heMaxp(thermo.he(Tmaxp, patchI));

                forAll(hep, faceI)
                {
                    hep[faceI] =
                        max(min(hep[faceI], heMaxp[faceI]), heMinp[faceI]);
                }
            }
        }
    }
    else  // T
    {
    scalarField& Tc = heT.ref();
        forAll(cells_, i)
        {
            label cellI = cells_[i];
            Tc[cellI] = max(min(Tc[cellI], Tmax_), Tmin_);
        }

        // handle boundaries in the case of 'all'
        if (selectionMode_ == smAll)
        {
            volScalarField::Boundary& bf = heT.boundaryFieldRef();

            forAll(bf, patchI)
            {
                fvPatchScalarField& Tp = bf[patchI];
                if (Tp.fixesValue())
                {
                    // not over-riding fixed conditions
                    continue;
                }

                forAll(Tp, faceI)
                {
                    Tp[faceI] = max(min(Tp[faceI], Tmax_), Tmin_);
                }
            }
        }
    }
}


void Foam::fv::temperatureLimitsConstraint::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fv::temperatureLimitsConstraint::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        coeffs_.readIfPresent("Tmin", Tmin_);
        coeffs_.readIfPresent("Tmax", Tmax_);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
