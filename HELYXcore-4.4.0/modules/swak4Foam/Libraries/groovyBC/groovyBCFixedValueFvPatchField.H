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

Class
    Foam::groovyBCFixedValueFvPatchField

Description
    Foam::groovyBCFixedValueFvPatchField

SourceFiles
    groovyBCFixedValueFvPatchField.C

Contributors/Copyright:
    2009-2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef groovyBCFixedValueFvPatchField_H
#define groovyBCFixedValueFvPatchField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchField.H"

#include "PatchValueExpressionDriver.H"

#include "groovyBCCommon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class groovyBCFixedValueFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class groovyBCFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>,
    public groovyBCCommon<Type>
{

protected:

    // Protected data

    PatchValueExpressionDriver driver_;

public:

    //- Runtime type information
    TypeName("groovyBCFixedValue");


    // Constructors

    //- Construct from patch and internal field
    groovyBCFixedValueFvPatchField
    (
        const fvPatch&,
        const DimensionedField<Type, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    groovyBCFixedValueFvPatchField
    (
        const fvPatch&,
        const DimensionedField<Type, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given groovyBCFixedValueFvPatchField onto a new patch
    groovyBCFixedValueFvPatchField
    (
        const groovyBCFixedValueFvPatchField<Type>&,
        const fvPatch&,
        const DimensionedField<Type, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Construct as copy
    groovyBCFixedValueFvPatchField
    (
        const groovyBCFixedValueFvPatchField<Type>&
    );

    //- Construct and return a clone
    virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
                (
                    new groovyBCFixedValueFvPatchField<Type>(*this)
                );
        }

    //- Construct as copy setting internal field reference
    groovyBCFixedValueFvPatchField
    (
        const groovyBCFixedValueFvPatchField<Type>&,
        const DimensionedField<Type, volMesh>&
    );

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchField<Type>> clone
    (
        const DimensionedField<Type, volMesh>& iF
    ) const
        {
            return tmp<fvPatchField<Type>>
                (
                    new groovyBCFixedValueFvPatchField<Type>(*this, iF)
                );
        }


    // Member functions

    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();

    //- Write
    virtual void write(Ostream&) const;
};


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #endif

// ************************************************************************* //
