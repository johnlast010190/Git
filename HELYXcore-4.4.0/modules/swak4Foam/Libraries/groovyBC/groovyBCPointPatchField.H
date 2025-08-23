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

Class
    Foam::groovyBCPointPatchField

Description

SourceFiles
    groovyBCPointPatchField.C

Contributors/Copyright:
    2010-2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef groovyBCPointPatchField_H
#define groovyBCPointPatchField_H

#include "include/swak.H"

    #ifdef FOAM_POINTPATCHFIELD_HAS_FIVE_TEMPLATE_PARAMETERS

#   include "MixedPointPatchField.H"
    // Predefine extra template paramters
#   define mixedPointPatchFieldType MixedPointPatchField<pointPatchField,pointMesh,pointPatch,DummyMatrix,Type>
#   define pointPatchFieldMapper PointPatchFieldMapper

#else

#ifdef FOAM_NO_MIXED_POINT_PATCH
#include "fields/pointPatchFields/basic/value/valuePointPatchField.H"
#    define mixedPointPatchFieldType valuePointPatchField<Type>
#else
#    include "mixedPointPatchField.H"
#    define mixedPointPatchFieldType mixedPointPatchField<Type>
#endif

#endif

#include "PatchValueExpressionDriver.H"

#include "groovyBCCommon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class groovyBCPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class groovyBCPointPatchField
:
    public mixedPointPatchFieldType,
    public groovyBCCommon<Type>
{

protected:

    // Protected data

    PatchValueExpressionDriver driver_;

    string nullValue();

public:

    //- Runtime type information
    TypeName("groovyBC");


    // Constructors

    //- Construct from patch and internal field
    groovyBCPointPatchField
    (
        const pointPatch&,
        const DimensionedField<Type, pointMesh>&
    );

    //- Construct from patch, internal field and dictionary
    groovyBCPointPatchField
    (
        const pointPatch&,
        const DimensionedField<Type, pointMesh>&,
        const dictionary&
    );

    //- Construct by mapping given patchField<Type> onto a new patch
    groovyBCPointPatchField
    (
        const groovyBCPointPatchField<Type>&,
        const pointPatch&,
        const DimensionedField<Type, pointMesh>&,
        const pointPatchFieldMapper&
    );

    //- Construct and return a clone
    virtual autoPtr<pointPatchField<Type>> clone() const
        {
            return autoPtr<pointPatchField<Type>>
                (
                    new groovyBCPointPatchField<Type>
                    (
                        *this
                    )
                );
        }

    //- Construct as copy setting internal field reference
    groovyBCPointPatchField
    (
        const groovyBCPointPatchField<Type>&,
        const DimensionedField<Type, pointMesh>&
    );

    //- Construct and return a clone setting internal field reference
    virtual autoPtr<pointPatchField<Type>> clone
    (
        const DimensionedField<Type, pointMesh>& iF
    ) const
        {
            return autoPtr<pointPatchField<Type>>
                (
                    new groovyBCPointPatchField<Type>
                    (
                        *this,
                        iF
                    )
                );
        }


    // Member functions

    // Evaluation functions

    //- Update the patch field
    virtual void updateCoeffs();


    //- Write
    virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
