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
    (c) 2011-2024 Engys Ltd.

Class
    Foam::blendingSource

Description
    Abstract base class for blending functions

SourceFiles
    blendingSource.C

\*---------------------------------------------------------------------------*/

#ifndef blendingSource_H
#define blendingSource_H

#include "fieldBlendingFactor/blendingFunctions/blendingFunction/blendingFunction.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                          Class blendingSource Declaration
\*---------------------------------------------------------------------------*/

class blendingSource
{
protected:

    // Protected Data

        //- Value to coefficient translator
        autoPtr<blendingFunction> ramp_;


private:

    // Private Data

        //- Output switch
        Switch writeBlendingField_;


    // Protected Member Functions

        //- Calculate and apply blending factor
        tmp<surfaceScalarField> convertValueToBlendingFactor
        (
            const surfaceScalarField& valueField,
            bool verbose = true
        );


protected:

    // Protected Data
        //- object registry
        const objectRegistry& obr_;

        //- const reference to mesh
        const fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("blendingSource");


    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        blendingSource,
        dictionary,
        (
            const objectRegistry& obr,
            const fvMesh& mesh,
            const dictionary& blendingSourceDict
        ),
        (obr, mesh, blendingSourceDict)
    );


    // Constructors

        //- Construct from components
        blendingSource
        (
            const objectRegistry& obr,
            const fvMesh& mesh,
            const dictionary& data
        );

        //- Disallow default bitwise copy construction
        blendingSource(const blendingSource&) = delete;


    // Selectors

        //- Return an autoPtr to the selected blendingSource
        static autoPtr<blendingSource> New
        (
            const objectRegistry& obr,
            const fvMesh& mesh,
            const dictionary& blendingSourceDict
        );


    //- Destructor
    virtual ~blendingSource()
    {}


    // Member Functions

        //- return blending field from this
        tmp<surfaceScalarField> blendingFactor();

        // write blending source field to disk
        void write();

        //- calculate source field
        virtual tmp<surfaceScalarField> sourceField() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const blendingSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
