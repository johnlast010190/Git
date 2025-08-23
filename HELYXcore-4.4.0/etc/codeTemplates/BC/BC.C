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
    (c) YEAR Engys Ltd

\*---------------------------------------------------------------------------*/

#include "CONSTRUCT.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::scalar Foam::CLASS::t() const
{
    return this->db().time().timeOutputValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::CLASS::
CONSTRUCT
(
    const fvPatch& p,
    const DimensionedField<TYPE, volMesh>& iF
)
:
    PARENT(p, iF),
    scalarData_(0.0),
    data_(Zero),
    fieldData_(p.size(), Zero),
    timeVsData_(),
    wordData_("wordDefault"),
    labelData_(-1),
    boolData_(false)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::CLASS::
CONSTRUCT
(
    const fvPatch& p,
    const DimensionedField<TYPE, volMesh>& iF,
    const dictionary& dict
)
:
    PARENT(p, iF),
    scalarData_(dict.lookup<scalar>("scalarData")),
    data_(pTraits<TYPE>(dict.lookup("data"))),
    fieldData_("fieldData", dict, p.size()),
    timeVsData_(Function1<TYPE>::New("timeVsData", dict)),
    wordData_(dict.lookupOrDefault<word>("wordName", "wordDefault")),
    labelData_(-1),
    boolData_(false)
{
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;

    this->refValue() = FIELD("fieldData", dict, p.size());
    FVPATCHF::operator=(this->refValue());

    PARENT::evaluate();

    /*
    //Initialise with the value entry if evaluation is not possible
    FVPATCHF::operator=
    (
        FIELD("value", dict, p.size())
    );
    this->refValue() = *this;
    */
}


template<class Type>
Foam::CLASS::
CONSTRUCT
(
    const CLASS& ptf,
    const fvPatch& p,
    const DimensionedField<TYPE, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    PARENT(ptf, p, iF, mapper),
    scalarData_(ptf.scalarData_),
    data_(ptf.data_),
    fieldData_(ptf.fieldData_, mapper),
    timeVsData_(ptf.timeVsData_, false),
    wordData_(ptf.wordData_),
    labelData_(-1),
    boolData_(ptf.boolData_)
{}


template<class Type>
Foam::CLASS::
CONSTRUCT
(
    const CLASS& ptf
)
:
    PARENT(ptf),
    scalarData_(ptf.scalarData_),
    data_(ptf.data_),
    fieldData_(ptf.fieldData_),
    timeVsData_(ptf.timeVsData_, false),
    wordData_(ptf.wordData_),
    labelData_(-1),
    boolData_(ptf.boolData_)
{}


template<class Type>
Foam::CLASS::
CONSTRUCT
(
    const CLASS& ptf,
    const DimensionedField<TYPE, volMesh>& iF
)
:
    PARENT(ptf, iF),
    scalarData_(ptf.scalarData_),
    data_(ptf.data_),
    fieldData_(ptf.fieldData_),
    timeVsData_(ptf.timeVsData_, false),
    wordData_(ptf.wordData_),
    labelData_(-1),
    boolData_(ptf.boolData_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::CLASS::autoMap
(
    const fvPatchFieldMapper& m
)
{
    PARENT::autoMap(m);
    fieldData_.autoMap(m);
}


template<class Type>
void Foam::CLASS::rmap
(
    const FVPATCHF& ptf,
    const labelList& addr
)
{
    PARENT::rmap(ptf, addr);

    const CLASS& tiptf =
        refCast<const CLASS>(ptf);

    fieldData_.rmap(tiptf.fieldData_, addr);
}


template<class Type>
void Foam::CLASS::reset
(
    const FVPATCHF& ptf
)
{
    PARENT::reset(ptf);

    const CLASS& tiptf =
        refCast<const CLASS>(ptf);

    fieldData_.reset(tiptf.fieldData_);
}


template<class Type>
void Foam::CLASS::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    PARENT::autoMapGIB(mapper);
    mapper.map(fieldData_, defaultValue);
}


template<class Type>
void Foam::CLASS::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    PARENT::forceAssign
    (
        data_
      + fieldData_
      + scalarData_*timeVsData_->value(t())
    );

    const scalarField& phip =
        this->patch().template lookupPatchField<surfaceScalarField, scalar>
        (
            "phi"
        );
    this->valueFraction() = 1.0 - pos(phip);

    PARENT::updateCoeffs();
}


template<class Type>
void Foam::CLASS::write
(
    Ostream& os
) const
{
    FVPATCHF::write(os);
    os.writeEntry("scalarData", scalarData_);
    os.writeEntry("data", data_):
    fieldData_.writeEntry("fieldData", os);
    timeVsData_->writeData(os);
    os.writeEntry("wordData", wordData_);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        FVPATCHF,
        CLASS
    );
}

// ************************************************************************* //
