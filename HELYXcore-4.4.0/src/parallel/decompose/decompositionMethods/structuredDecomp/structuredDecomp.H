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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::structuredDecomp

Description
    Decomposition by walking out decomposition of patch cells mesh.

SourceFiles
    structuredDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef structuredDecomp_H
#define structuredDecomp_H

#include "decompositionMethod/decompositionMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class structuredDecomp Declaration
\*---------------------------------------------------------------------------*/

class structuredDecomp
:
    public decompositionMethod
{
    // Private data

        dictionary methodDict_;

        wordReList patches_;

        autoPtr<decompositionMethod> method_;


public:

    //- Runtime type information
    TypeName("structured");


    // Constructors

        //- Construct given the decomposition dictionary
        structuredDecomp(const dictionary& decompositionDict);

        //- Disallow default bitwise assignment
        structuredDecomp(const structuredDecomp&) = delete;


    //- Destructor
    virtual ~structuredDecomp()
    {}


    // Member Functions

        //- Is method parallel aware (i.e. does it synchronize domains across
        //  proc boundaries)
        virtual bool parallelAware() const;

        //- Return for every coordinate the wanted processor number. Use the
        //  mesh connectivity (if needed)
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& points,
            const scalarField& pointWeights
        );

        //- Return for every coordinate the wanted processor number. Explicitly
        //  provided connectivity - does not use mesh_.
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights
        );


    // Member Operators

        //- Disallow default bitwise copy construction
        void operator=(const structuredDecomp&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
