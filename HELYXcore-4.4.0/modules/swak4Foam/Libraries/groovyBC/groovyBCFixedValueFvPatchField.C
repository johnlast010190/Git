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

Contributors/Copyright:
    2009-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "groovyBCFixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
groovyBCFixedValueFvPatchField<Type>::groovyBCFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    groovyBCCommon<Type>(false),
    driver_(this->patch())
{
    if (debug) {
        Info<< "groovyBCFixedValueFvPatchField<Type>::groovyBCFixedValueFvPatchField 1" << endl;
    }
}


template<class Type>
groovyBCFixedValueFvPatchField<Type>::groovyBCFixedValueFvPatchField
(
    const groovyBCFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    groovyBCCommon<Type>(ptf),
    driver_(this->patch(),ptf.driver_)
{
    if (debug) {
        Info<< "groovyBCFixedValueFvPatchField<Type>::groovyBCFixedValueFvPatchField 2" << endl;
    }
}


template<class Type>
groovyBCFixedValueFvPatchField<Type>::groovyBCFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    groovyBCCommon<Type>(dict,false),
    driver_(dict,this->patch(),this->db())
{
    if (debug) {
        Info<< "groovyBCFixedValueFvPatchField<Type>::groovyBCFixedValueFvPatchField 3" << endl;
    }

    driver_.readVariablesAndTables(dict);

    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        this->forceAssign(this->patchInternalField());

        WarningIn(
            "groovyBCFixedValueFvPatchField<Type>::groovyBCFixedValueFvPatchField"
            "("
            "const fvPatch& p,"
            "const DimensionedField<Type, volMesh>& iF,"
            "const dictionary& dict"
            ")"
        ) << "No value defined for "
#ifdef FOAM_NO_DIMENSIONEDINTERNAL_IN_GEOMETRIC
            << this->internalField().name()
#else
            << this->dimensionedInternalField().name()
#endif
            << " on " << this->patch().name() << " therefore would be undefined "
            << "and set to the internal field next to the patch"
            << endl;
    }

    if (this->evaluateDuringConstruction()) {
        // make sure that this works with potentialFoam or other solvers that don't evaluate the BCs
        this->evaluate();
    } else {
        // original does nothing
    }
}


template<class Type>
groovyBCFixedValueFvPatchField<Type>::groovyBCFixedValueFvPatchField
(
    const groovyBCFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    groovyBCCommon<Type>(ptf),
    driver_(this->patch(),ptf.driver_)
{
    if (debug) {
        Info<< "groovyBCFixedValueFvPatchField<Type>::groovyBCFixedValueFvPatchField 4" << endl;
    }
}


template<class Type>
groovyBCFixedValueFvPatchField<Type>::groovyBCFixedValueFvPatchField
(
    const groovyBCFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    groovyBCCommon<Type>(ptf),
    driver_(this->patch(),ptf.driver_)
{
    if (debug) {
        Info<< "groovyBCFixedValueFvPatchField<Type>::groovyBCFixedValueFvPatchField 5" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void groovyBCFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (debug) {
        Info<< "groovyBCFixedValueFvPatchField<Type>::updateCoeffs" << endl;
        Info<< "Value: " << this->valueExpression_ << endl;
        Info<< "Gradient: " << this->gradientExpression_ << endl;
        Info<< "Fraction: " << this->fractionExpression_ << endl;
        Info<< "Variables: ";
        driver_.writeVariableStrings(Info) << endl;
    }
    if (this->updated())
    {
        return;
    }

    if (debug) {
        Info<< "groovyBCFixedValueFvPatchField<Type>::updateCoeffs - updating" << endl;
    }

    driver_.clearVariables();

    this->forceAssign(driver_.evaluate<Type>(this->valueExpression_));

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void groovyBCFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    if (debug) {
        Info<< "groovyBCFixedValueFvPatchField<Type>::write" << endl;
    }
    fixedValueFvPatchField<Type>::write(os);
    groovyBCCommon<Type>::write(os);

    driver_.writeCommon(os,this->debug_ || debug);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
