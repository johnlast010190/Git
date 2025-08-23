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
    (c) 2024 Engys Ltd.

Class
    Foam::lcsSpeciesSourcePluginFunction

Description

SourceFiles
    lcsSpeciesSourcePluginFunction.C

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef lcsSpeciesSourcePluginFunction_H
#define lcsSpeciesSourcePluginFunction_H

#include "LagrangianCloudSourcePluginFunction.H"
#include "include/swakCloudTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class lcsSpeciesSourcePluginFunction Declaration
\*---------------------------------------------------------------------------*/

class lcsSpeciesSourcePluginFunction
:
    public LagrangianCloudSourcePluginFunction
{
protected:

    // Protected Data

        word speciesName_;

        typedef DimensionedField<scalar,volMesh> dimScalarField;


    // Protected Member Functions

        virtual autoPtr<dimScalarField> internalEvaluate(const label index);

        virtual label getIndex(wordList& specList);


public:

    //- Runtime type information
    TypeName("lcsSpeciesSourcePluginFunction");


    // Constructors

        //- Construct from components
        lcsSpeciesSourcePluginFunction
        (
            const FieldValueExpressionDriver& parentDriver,
            const word& name
        );

        //- Disallow default bitwise copy construction
        lcsSpeciesSourcePluginFunction
        (
            const lcsSpeciesSourcePluginFunction&
        ) = delete;


    //- Destructor
    virtual ~lcsSpeciesSourcePluginFunction() {}


    // Member Functions

        //- Perform the evaluation
        void doEvaluation();

        void setArgument(label index,const word &);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const lcsSpeciesSourcePluginFunction&) = delete;
};

#ifdef FOAM_REACTINGCLOUD_TEMPLATED
#define getSpeciesIndex(neededType,usedType)                                \
    if (hasCloudAs<neededType>()) {                                      \
        const speciesTable &spec=                                       \
            getCloudAs<neededType,usedType>().mcCarrierThermo ().species(); \
        speciesList=spec;                                               \
        if (spec.found(speciesName_)) {                                 \
            return spec[speciesName_];                                  \
        }                                                               \
    }
#else
#define getSpeciesIndex(neededType,usedType)                            \
    if (hasCloudAs<neededType>()) {                                      \
        const speciesTable &spec=                                       \
            getCloudAs<neededType,usedType>().thermo().carrier().species(); \
        speciesList=spec;                                               \
        if (spec.found(speciesName_)) {                                 \
            return spec[speciesName_];                                  \
        }                                                               \
    }
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
