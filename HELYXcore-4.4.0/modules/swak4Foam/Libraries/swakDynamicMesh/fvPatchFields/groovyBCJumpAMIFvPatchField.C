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
    2011, 2013-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "include/swak.H"

#if defined(FOAM_HAS_AMI_INTERFACE)

#include "groovyBCJumpAMIFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
groovyBCJumpAMIFvPatchField<Type>::groovyBCJumpAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    jumpCyclicAMIFvPatchField<Type>(p, iF),
    driver_(this->patch()),
    jumpExpression_("0")
{
    if (debug) {
        Info<< "groovyBCJumpAMIFvPatchField<Type>::groovyBCJumpAMIFvPatchField 1" << endl;
    }
}


template<class Type>
groovyBCJumpAMIFvPatchField<Type>::groovyBCJumpAMIFvPatchField
(
    const groovyBCJumpAMIFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    jumpCyclicAMIFvPatchField<Type>(ptf, p, iF, mapper),
    driver_(this->patch(),ptf.driver_),
    jumpExpression_(ptf.jumpExpression_)
{
    if (debug) {
        Info<< "groovyBCJumpAMIFvPatchField<Type>::groovyBCJumpAMIFvPatchField 2" << endl;
    }
}


template<class Type>
groovyBCJumpAMIFvPatchField<Type>::groovyBCJumpAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    jumpCyclicAMIFvPatchField<Type>(p, iF),
    driver_(dict,this->patch(), this->db()),
    jumpExpression_(
        dict.lookup("jumpExpression"),
        dict
    )
{
    if (debug) {
        Info<< "groovyBCJumpAMIFvPatchField<Type>::groovyBCJumpAMIFvPatchField 3" << endl;
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
        this->evaluate(
#ifdef FOAM_PSTREAM_COMMSTYPE_IS_ENUMCLASS
            Pstream::commsTypes::blocking
#else
            Pstream::blocking
#endif
        );
    }
}


template<class Type>
groovyBCJumpAMIFvPatchField<Type>::groovyBCJumpAMIFvPatchField
(
    const groovyBCJumpAMIFvPatchField<Type>& ptf
)
:
    jumpCyclicAMIFvPatchField<Type>(ptf),
    driver_(this->patch(),ptf.driver_),
    jumpExpression_(ptf.jumpExpression_)
{
    if (debug) {
        Info<< "groovyBCJumpAMIFvPatchField<Type>::groovyBCJumpAMIFvPatchField 4" << endl;
    }
}


template<class Type>
groovyBCJumpAMIFvPatchField<Type>::groovyBCJumpAMIFvPatchField
(
    const groovyBCJumpAMIFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    jumpCyclicAMIFvPatchField<Type>(ptf, iF),
    driver_(this->patch(),ptf.driver_),
    jumpExpression_(ptf.jumpExpression_)
{
    if (debug) {
        Info<< "groovyBCJumpAMIFvPatchField<Type>::groovyBCJumpAMIFvPatchField 5" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type>> groovyBCJumpAMIFvPatchField<Type>::jump() const
{
    if (debug) {
        Info<< "groovyBCJumpAMIFvPatchField<Type>::jump() with "
            << jumpExpression_ << endl;
    }

    if (this->cyclicAMIPatch().owner()) {
        PatchValueExpressionDriver &driver
            = const_cast<PatchValueExpressionDriver &>(driver_);

        driver.clearVariables();

        tmp<Field<Type>> tjf(
            new Field<Type>(this->size())
        );
        Field<Type> &jf=const_cast<Field<Type>&>(tjf());

        jf = driver.evaluate<Type>(this->jumpExpression_);

        return jf;
    } else {
       const groovyBCJumpAMIFvPatchField& nbrPatch =
           refCast<const groovyBCJumpAMIFvPatchField<Type>>(
               this->nbrPatchField()
           );
       if (this->cyclicAMIPatch().applyLowWeightCorrection()) {
           return this->cyclicAMIPatch().interpolate(
               nbrPatch.jump(),
               Field<Type>(this->size(), Zero)
           );
       } else {
           return this->cyclicAMIPatch().interpolate(nbrPatch.jump());
       }
    }
}

template<class Type>
void groovyBCJumpAMIFvPatchField<Type>::write(Ostream& os) const
{
    if (debug) {
        Info<< "groovyBCJumpAMIFvPatchField<Type>::write" << endl;
    }
    fvPatchField<Type>::write(os);
    if (!this->overridesConstraint())
    {
        os.writeEntry("patchType", "cyclic");
    }
    if (this->cyclicAMIPatch().owner()) {
        os.writeEntry("jumpValue", jump());
    }
    os.writeEntry("jumpExpression", jumpExpression_);
    driver_.writeCommon(os,debug);

    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
