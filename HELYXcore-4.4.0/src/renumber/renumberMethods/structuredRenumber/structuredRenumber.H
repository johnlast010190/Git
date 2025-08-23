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
    (c) 2012-2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::structuredRenumber

Description
    Renumbering according to mesh layers.
    depthFirst = true:
        first column gets ids 0..nLayer-1,
        second nLayers..2*nLayers-1 etc.
    depthFirst = false:
        first layer gets ids 0,1,2 etc.

SourceFiles
    structuredRenumber.C

\*---------------------------------------------------------------------------*/

#ifndef structuredRenumber_H
#define structuredRenumber_H

#include "renumberMethod/renumberMethod.H"
#include "meshStructure/topoDistanceData.H"
#include "primitives/bools/Switch/Switch.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class structuredRenumber Declaration
\*---------------------------------------------------------------------------*/

class structuredRenumber
:
    public renumberMethod
{
public:

    // Public classes

        //- Less function class that can be used for sorting according to
        //  column and layer
        class layerLess
        {
            const Switch depthFirst_;
            const labelList& order_;
            const List<topoDistanceData>& distance_;

        public:

            layerLess
            (
                const Switch depthFirst,
                const labelList& order,
                const List<topoDistanceData>& distance
            )
            :
                depthFirst_(depthFirst),
                order_(order),
                distance_(distance)
            {}

            bool operator()(const label a, const label b);
        };


    // Private data

        const dictionary methodDict_;

        const wordReList patches_;

        const label nLayers_;

        const Switch depthFirst_;

        const autoPtr<renumberMethod> method_;

        const Switch reverse_;


public:

    //- Runtime type information
    TypeName("structured");


    // Constructors

        //- Construct given the renumber dictionary
        structuredRenumber(const dictionary& renumberDict);

        //- Disallow default bitwise assignment
        structuredRenumber(const structuredRenumber&) = delete;


    //- Destructor
    virtual ~structuredRenumber()
    {}


    // Member Functions

        //- Return the order in which cells need to be visited, i.e.
        //  from ordered back to original cell label.
        //  This is only defined for geometric renumberMethods.
        virtual labelList renumber(const pointField&) const
        {
            NotImplemented;
        }

        //- Return the order in which cells need to be visited, i.e.
        //  from ordered back to original cell label.
        //  Use the mesh connectivity (if needed)
        virtual labelList renumber
        (
            const polyMesh& mesh,
            const pointField& cc
        ) const;

        //- Return the order in which cells need to be visited, i.e.
        //  from ordered back to original cell label.
        //  The connectivity is equal to mesh.cellCells() except
        //  - the connections are across coupled patches
        virtual labelList renumber
        (
            const labelListList& cellCells,
            const pointField& cc
        ) const
        {
            NotImplemented;
        }


    // Member Operators

        //- Disallow default bitwise copy construction
        void operator=(const structuredRenumber&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
