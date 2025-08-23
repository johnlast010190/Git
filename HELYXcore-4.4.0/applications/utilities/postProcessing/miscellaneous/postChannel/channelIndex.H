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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd

Class
    Foam::channelIndex

Description
    Does averaging of fields over layers of cells. Assumes layered mesh.

SourceFiles
    channelIndex.C

\*---------------------------------------------------------------------------*/

#ifndef channelIndex_H
#define channelIndex_H

#include "regionSplit/regionSplit.H"
#include "primitives/direction/direction.H"
#include "fields/Fields/scalarField/scalarField.H"
#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class channelIndex Declaration
\*---------------------------------------------------------------------------*/

class channelIndex
{

    // Private data

        static const NamedEnum<vector::components, 3> vectorComponentsNames_;

        //- Is mesh symmetric
        const bool symmetric_;

        //- Direction to sort
        const direction dir_;

        //- Per cell the global region
        autoPtr<regionSplit> cellRegion_;

        //- Per global region the number of cells (scalarField so we can use
        //  field algebra)
        scalarField regionCount_;

        //- From sorted region back to unsorted global region
        labelList sortMap_;

        //- Sorted component of cell centres
        scalarField y_;



    // Private Member Functions

        void walkOppositeFaces
        (
            const polyMesh& mesh,
            const labelList& startFaces,
            boolList& blockedFace
        );

        void calcLayeredRegions
        (
            const polyMesh& mesh,
            const labelList& startFaces
        );


public:

    // Constructors

        //- Construct from dictionary
        channelIndex(const polyMesh&, const dictionary&);

        //- Construct from supplied starting faces
        channelIndex
        (
            const polyMesh& mesh,
            const labelList& startFaces,
            const bool symmetric,
            const direction dir
        );

        //- Disallow default bitwise copy construction
        channelIndex(const channelIndex&) = delete;


    // Member Functions

        // Access

            //- Sum field per region
            template<class T>
            Field<T> regionSum(const Field<T>& cellField) const;

            //- Collapse a field to a line
            template<class T>
            Field<T> collapse
            (
                const Field<T>& vsf,
                const bool asymmetric=false
            ) const;

            //- Return the field of Y locations from the cell centres
            const scalarField& y() const
            {
                return y_;
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const channelIndex&) = delete;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "channelIndexTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
