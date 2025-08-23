/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) ICE Stroemungsfoschungs GmbH
    (c) 1991-2010 OpenCFD Ltd.

Contributors/Copyright:
    2011, 2013-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "groovyFixedNormalSlipPointPatchField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

#include "groovyBCCommon.H"

#ifdef FOAM_POINTPATCHFIELD_HAS_FIVE_TEMPLATE_PARAMETERS
#include "PointPatchFieldMapper.H"
#else
#include "fields/pointPatchFields/pointPatchField/pointPatchFieldMapper.H"
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::groovyFixedNormalSlipPointPatchField<Type>::groovyFixedNormalSlipPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    SlipPointPatchFieldType(p, iF),
    fixedValueExpression_(""),
    normalExpression_("toPoint(normal())"),
    driver_(groovyBCCommon<Type>::getFvPatch(this->patch()))
{}


template<class Type>
Foam::groovyFixedNormalSlipPointPatchField<Type>::groovyFixedNormalSlipPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    SlipPointPatchFieldType(p, iF),
    fixedValueExpression_(
    dict.lookup("fixedValueExpression"),
    dict
    ),
    normalExpression_(
    dict.lookup("normalExpression"),
    dict
    ),
    driver_(dict,groovyBCCommon<Type>::getFvPatch(this->patch()),this->db())
{
}


template<class Type>
Foam::groovyFixedNormalSlipPointPatchField<Type>::groovyFixedNormalSlipPointPatchField
(
    const groovyFixedNormalSlipPointPatchField& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    SlipPointPatchFieldType(ptf, p, iF, mapper),
    fixedValueExpression_(ptf.fixedValueExpression_),
    normalExpression_(ptf.normalExpression_),
    driver_(groovyBCCommon<Type>::getFvPatch(this->patch()),ptf.driver_)
{}


template<class Type>
Foam::groovyFixedNormalSlipPointPatchField<Type>::groovyFixedNormalSlipPointPatchField
(
    const groovyFixedNormalSlipPointPatchField& tppsf
)
:
    SlipPointPatchFieldType(tppsf),
    fixedValueExpression_(tppsf.fixedValueExpression_),
    normalExpression_(tppsf.normalExpression_),
    driver_(groovyBCCommon<Type>::getFvPatch(this->patch()),tppsf.driver_)
{}


template<class Type>
Foam::groovyFixedNormalSlipPointPatchField<Type>::groovyFixedNormalSlipPointPatchField
(
    const groovyFixedNormalSlipPointPatchField& tppsf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    SlipPointPatchFieldType(tppsf, iF),
    fixedValueExpression_(tppsf.fixedValueExpression_),
    normalExpression_(tppsf.normalExpression_),
    driver_(groovyBCCommon<Type>::getFvPatch(this->patch()),tppsf.driver_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::groovyFixedNormalSlipPointPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    driver_.clearVariables();

    Field<vector> n(driver_.evaluate<vector>(this->normalExpression_,true));
    n/=mag(n); // normalize
    Field<Type> val(driver_.evaluate<Type>(this->fixedValueExpression_,true));

    tmp<Field<Type>> tvalues=transform(I - n*n, this->patchInternalField())+transform(n*n,val);

    // Get internal field to insert values into
    Field<Type>& iF =
#ifdef FOAM_NO_DIMENSIONEDINTERNAL_IN_GEOMETRIC
        const_cast<Field<Type>&>(this->internalField().field());
#else
    const_cast<Field<Type>&>(this->internalField());
#endif

    this->setInternalField(iF, tvalues());
}

template<class Type>
void Foam::groovyFixedNormalSlipPointPatchField<Type>::write(Ostream& os) const
{
    SlipPointPatchFieldType::write(os);

    os.writeEntry("fixedValueExpression", fixedValueExpression_);
    os.writeEntry("normalExpression", normalExpression_);
    os.writeEntry("value", this->patchInternalField());

    if (this->fixedValueExpression_!="") {
        os.writeEntry
        (
            "fixedValue",
            const_cast<PatchValueExpressionDriver&>(driver_).evaluate<Type>
            (
                this->fixedValueExpression_,
                true
            )
        );
    }
    os.writeEntry
    (
        "normal",
        const_cast<PatchValueExpressionDriver&>(driver_).evaluate<Type>
        (
            this->normalExpression_,
            true
        )
    );

    driver_.writeCommon(os,debug);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
