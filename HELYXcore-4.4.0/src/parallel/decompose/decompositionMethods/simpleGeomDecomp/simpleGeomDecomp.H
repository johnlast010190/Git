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
    (c) 2011 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::simpleGeomDecomp

Description

SourceFiles
    simpleGeomDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef simpleGeomDecomp_H
#define simpleGeomDecomp_H

#include "geomDecomp/geomDecomp.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class simpleGeomDecomp Declaration
\*---------------------------------------------------------------------------*/

class simpleGeomDecomp
:
    public geomDecomp
{
    // Private Member Functions

        void assignToProcessorGroup(labelList&, const label) const;

        void assignToProcessorGroup
        (
            labelList& processorGroup,
            const label nProcGroup,
            const labelList& indices,
            const scalarField& weights,
            const scalar summedWeights
        ) const;

        labelList decomposeOneProc(const pointField& points) const;

        labelList decomposeOneProc
        (
            const pointField& points,
            const scalarField& weights
        ) const;


public:

    //- Runtime type information
    TypeName("simple");


    // Constructors

        //- Construct given the decomposition dictionary
        simpleGeomDecomp(const dictionary& decompositionDict);

        //- Disallow default bitwise assignment
        simpleGeomDecomp(const simpleGeomDecomp&) = delete;


    //- Destructor
    virtual ~simpleGeomDecomp()
    {}


    // Member Functions

        virtual bool parallelAware() const
        {
            // simpleDecomp sends all points to the master which does
            // the decomposition.
            return true;
        }

        virtual labelList decompose(const pointField&);

        virtual labelList decompose(const pointField&, const scalarField&);

        virtual labelList decompose(const polyMesh&, const pointField& points)
        {
            return decompose(points);
        }

        virtual labelList decompose
        (
            const polyMesh&,
            const pointField& points,
            const scalarField& pointWeights
        )
        {
            return decompose(points, pointWeights);
        }

        //- Explicitly provided connectivity
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights
        )
        {
            return decompose(cc, cWeights);
        }


    // Member Operators

        //- Disallow default bitwise copy construction
        void operator=(const simpleGeomDecomp&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
