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
    (c) 2012-2023 OpenFOAM Foundation
    (c) 2017-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fixedTemperatureConstraint.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(fixedTemperatureConstraint, 0);
        addToRunTimeSelectionTable
        (
            option,
            fixedTemperatureConstraint,
            dictionary
        );
    }

    template<>
    const char* NamedEnum<fv::fixedTemperatureConstraint::temperatureMode, 3>::
    names[] =
    {
        "uniform",
        "lookup",
        "spatial"
    };
}

const Foam::NamedEnum<Foam::fv::fixedTemperatureConstraint::temperatureMode, 3>
    Foam::fv::fixedTemperatureConstraint::temperatureModeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::fixedTemperatureConstraint::fixedTemperatureConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    coorFramePtr_(nullptr)
{
    fixedTemperatureConstraint::read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::fixedTemperatureConstraint::sourceFields(wordList& fieldNames)
{
    TName_ = coeffs_.lookupOrDefault<word>("T", "T");
    thermoPtr_ =
        obr_.lookupObjectRefPtr<basicThermo>(basicThermo::dictName);

    fieldNames.setSize
    (
        1, thermoPtr_ ? thermoPtr_->heT().name() : TName_
    );
}


void Foam::fv::fixedTemperatureConstraint::constrain
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    const scalar t = mesh_.time().value();

    scalarField Thevalue(cells_.size(), 0.0);
    switch (mode_)
    {
        case tmUniform:
        {
            Thevalue = scalarField(cells_.size(), Tfun1_().value(t));
            break;
        }
        case tmLookup:
        {
            const volScalarField& T =
                mesh_.lookupObject<volScalarField>(TName_);
            Thevalue = scalarField(T, cells_);
            break;
        }
        case tmSpatial:
        {
            vectorField C(mesh_.C().internalField(), cells_);
            if (coorFramePtr_)
            {
                C = coorFramePtr_->coorSys().localPosition(C);
            }
            Thevalue = Tfun1_().value(C);
            break;
        }
    }

    if (thermoPtr_ && thermoPtr_->calculatesTFromhe())
    {
        const scalarField pCells(thermoPtr_->p(), cells_);
        Thevalue = thermoPtr_->he(pCells, Thevalue, cells_);
    }

    // Optionally apply only a fraction of the constraint
    if (fraction_.valid())
    {
        eqn.setValues
        (
            cells_,
            Thevalue,
            scalarList(cells_.size(), fraction_->value(t))
        );
    }
    else
    {
        eqn.setValues(cells_, Thevalue);
    }
}


bool Foam::fv::fixedTemperatureConstraint::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        mode_ = temperatureModeNames_.read(coeffs_.lookup("mode"));

        if (mode_ != tmLookup)
        {
            Tfun1_ = Function1<scalar>::New("temperature", coeffs_);
        }

        if (coeffs_.found("referenceFrame"))
        {
            coorFramePtr_ = coordinateFrame::lookupNew(mesh_, coeffs_);
        }

        fraction_ =
            coeffs_.found("fraction")
          ? Function1<scalar>::New("fraction", coeffs_)
          : autoPtr<Function1<scalar>>();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
