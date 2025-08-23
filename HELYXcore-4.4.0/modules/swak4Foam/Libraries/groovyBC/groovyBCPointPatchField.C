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
    (c) 1991-2008 OpenCFD Ltd.

Contributors/Copyright:
    2010-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "groovyBCPointPatchField.H"

#ifdef FOAM_POINTPATCHFIELD_HAS_FIVE_TEMPLATE_PARAMETERS
#include "PointPatchFieldMapper.H"
#else
#include "fields/pointPatchFields/pointPatchField/pointPatchFieldMapper.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
groovyBCPointPatchField<Type>::groovyBCPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    mixedPointPatchFieldType(p, iF),
    groovyBCCommon<Type>(false,true),
    driver_(groovyBCCommon<Type>::getFvPatch(this->patch()))
{
#ifndef FOAM_NO_MIXED_POINT_PATCH
    this->refValue() = pTraits<Type>::zero;
    this->valueFraction() = 0.0;
#endif
}


template<class Type>
groovyBCPointPatchField<Type>::groovyBCPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    mixedPointPatchFieldType(p, iF),
    groovyBCCommon<Type>(dict,false,true),
    driver_(dict,groovyBCCommon<Type>::getFvPatch(this->patch()),this->db())
{
    driver_.readVariablesAndTables(dict);

#ifndef FOAM_NO_MIXED_POINT_PATCH
    if (dict.found("refValue")) {
        this->refValue() = Field<Type>("refValue", dict, p.size());
    } else {
        this->refValue() = this->patchInternalField();
    }
#endif

    if (dict.found("value"))
    {
        Field<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
#ifndef FOAM_NO_MIXED_POINT_PATCH
        Field<Type>::operator=(this->refValue());
#endif

        WarningIn(
            "groovyBCPointPatchField<Type>::groovyBCPointPatchField"
            "("
            "const pointPatch& p,"
            "const DimensionedField<Type, pointMesh>& iF,"
            "const dictionary& dict"
            ")"
        ) << "No value defined for "
#ifdef FOAM_NO_DIMENSIONEDINTERNAL_IN_GEOMETRIC
            << this->internalField().name()
#else
            << this->dimensionedInternalField().name()
#endif
            << " on " << this->patch().name() << " therefore using "
#ifndef FOAM_NO_MIXED_POINT_PATCH
            << "the internal field next to the patch"
#endif
            << endl;
    }

    //    this->refGrad() = pTraits<Type>::zero;

#ifndef FOAM_NO_MIXED_POINT_PATCH
    if (dict.found("valueFraction")) {
        this->valueFraction() = Field<scalar>("valueFraction", dict, p.size());
    } else {
        this->valueFraction() = 1;
    }
#endif

    if (this->evaluateDuringConstruction()) {
        // make sure that this works with potentialFoam or other solvers that don't evaluate the BCs
        this->evaluate();
    } else {
        // mixed-BC DOES NOT call evaluate during construction
    }
}


template<class Type>
groovyBCPointPatchField<Type>::groovyBCPointPatchField
(
    const groovyBCPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    mixedPointPatchFieldType
    (
        ptf,
        p,
        iF,
        mapper
    ),
    groovyBCCommon<Type>(ptf),
    driver_(groovyBCCommon<Type>::getFvPatch(this->patch()),ptf.driver_)
{
}


template<class Type>
groovyBCPointPatchField<Type>::groovyBCPointPatchField
(
    const groovyBCPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    mixedPointPatchFieldType(ptf, iF),
    groovyBCCommon<Type>(ptf),
    driver_(groovyBCCommon<Type>::getFvPatch(this->patch()),ptf.driver_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// Evaluate patch field
template<class Type>
void groovyBCPointPatchField<Type>::updateCoeffs()
{
    if (debug)
    {
        Info<< "groovyBCFvPatchField<Type>::updateCoeffs" << endl;
        Info<< "Value: " << this->valueExpression_ << endl;
        //        Info<< "Gradient: " << gradientExpression_ << endl;
        Info<< "Fraction: " << this->fractionExpression_ << endl;
        Info<< "Variables: ";
        driver_.writeVariableStrings(Info)  << endl;
    }
//     if (this->updated())
//     {
//         return;
//     }

    driver_.clearVariables();

#ifndef FOAM_NO_MIXED_POINT_PATCH
    this->refValue() = driver_.evaluate<Type>(this->valueExpression_,true);
    this->valueFraction() = driver_.evaluate<scalar>(this->fractionExpression_,true);
    //    this->refGrad() = driver_.evaluate<Type>(gradientExpression_,true);
#else
    Field<Type>::operator=
        (
            driver_.evaluate<Type>(this->valueExpression_,true)
        );
#endif

   mixedPointPatchFieldType::updateCoeffs();
}


// Write
template<class Type>
void groovyBCPointPatchField<Type>::write(Ostream& os) const
{
    mixedPointPatchFieldType::write(os);
    groovyBCCommon<Type>::write(os);

    this->writeEntry("value", os);

    driver_.writeCommon(os,this->debug_ || debug);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
