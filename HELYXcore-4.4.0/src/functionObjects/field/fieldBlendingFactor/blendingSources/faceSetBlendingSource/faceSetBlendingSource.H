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
    Foam::faceSetBlendingSource

Description
    faceSetBlendingSource blending source

SourceFiles
    faceSetBlendingSource.C

\*---------------------------------------------------------------------------*/

#ifndef faceSetBlendingSource_H
#define faceSetBlendingSource_H

#include "fieldBlendingFactor/blendingSources/blendingSource/blendingSource.H"
#include "sets/topoSets/faceSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingSources
{

/*---------------------------------------------------------------------------*\
                          Class faceSetBlendingSource Declaration
\*---------------------------------------------------------------------------*/

class faceSetBlendingSource
:
    public blendingSource
{
private:

    // Private Data

        //- set of faces to blend
        faceSet blendingFaceSet_;


public:

    //- Runtime type information
    TypeName("faceSet");


    // Constructors

        //- Construct from components
        faceSetBlendingSource
        (
            const objectRegistry& obr,
            const fvMesh& mesh,
            const dictionary& data
        );

        //- Disallow default bitwise copy construction
        faceSetBlendingSource(const faceSetBlendingSource&) = delete;


    //- Destructor
    virtual ~faceSetBlendingSource(){}


    // Member Functions

        //- source field calculation
        virtual tmp<surfaceScalarField> sourceField();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const faceSetBlendingSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blendingSources
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
