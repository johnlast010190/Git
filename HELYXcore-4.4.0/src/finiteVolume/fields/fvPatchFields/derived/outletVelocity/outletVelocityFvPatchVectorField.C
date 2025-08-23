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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/outletVelocity/outletVelocityFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::outletVelocityFvPatchVectorField::
outletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    relax_(0.2)
{}


Foam::outletVelocityFvPatchVectorField::
outletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    relax_(dict.lookupOrDefault<scalar>("relax", 0.2))
{
    if (dict.found("value"))
    {
        forceAssign(vectorField("value", dict, p.size()));
    }
    else
    {
        this->updateCoeffs();
    }
}


Foam::outletVelocityFvPatchVectorField::
outletVelocityFvPatchVectorField
(
    const outletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    relax_(ptf.relax_)
{}


Foam::outletVelocityFvPatchVectorField::
outletVelocityFvPatchVectorField
(
    const outletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    relax_(ptf.relax_)
{}


Foam::outletVelocityFvPatchVectorField::
outletVelocityFvPatchVectorField
(
    const outletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    relax_(ptf.relax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::outletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField::autoMap(m);
}


void Foam::outletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField::rmap(ptf, addr);
}


void Foam::outletVelocityFvPatchVectorField::reset
(
    const fvPatchVectorField& ptf
)
{
    fixedValueFvPatchVectorField::reset(ptf);
}


void Foam::outletVelocityFvPatchVectorField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchField::autoMapGIB(mapper);
}


Foam::tmp<Foam::vectorField>
Foam::outletVelocityFvPatchVectorField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    //- Assuming this BC is specified on outlet (phi > 0)
    tmp<vectorField> vict(new vectorField(this->size(), pTraits<vector>::one));

    return vict;
}


Foam::tmp<Foam::vectorField>
Foam::outletVelocityFvPatchVectorField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    //- Assuming this BC is specified on outlet (phi > 0)
    tmp<vectorField> vbct(new vectorField(this->size(), Zero));

    return vbct;
}


Foam::tmp<Foam::vectorField>
Foam::outletVelocityFvPatchVectorField::valueDivInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<vectorField>(new vectorField(this->size(), Zero));
}


Foam::tmp<Foam::vectorField>
Foam::outletVelocityFvPatchVectorField::valueDivBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<vectorField>(new vectorField(*this));
}


void Foam::outletVelocityFvPatchVectorField::updateCoeffs()
{
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::outletVelocityFvPatchVectorField::boundaryRelaxMatrix
(
    fvBlockMatrix<vector>& bEq
) const
{
    bEq.boundaryRelax(relax_, this->patch().index());
}


void Foam::outletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    this->writeEntry("value", os);
    writeEntryIfDifferent<scalar>(os, "relax", 0.2, relax_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        outletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
