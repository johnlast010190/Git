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
    (c) 2021-2024 Engys Ltd.

Class
    Foam::hexRef8DataList

Description
    Various for reading/decomposing/reconstructing/distributing refinement
    data for a list of meshes.

SourceFiles
    hexRef8DataList.C

\*---------------------------------------------------------------------------*/

#ifndef hexRef8DataList_H
#define hexRef8DataList_H

#include "primitives/ints/lists/labelIOList.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "containers/Lists/UPtrList/UPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyTopoChangeMap;
class polyDistributionMap;
class refinementHistory;
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class hexRef8DataList Declaration
\*---------------------------------------------------------------------------*/

class hexRef8DataList
{

private:

    // Private data

        //- mesh addressing (which proc has which mesh)
        labelList meshToProc_;

        labelListList meshNoInProc_;

        //- hexRef8 data
        PtrList<labelIOList> cellLevelPtrs_;

        PtrList<labelIOList> pointLevelPtrs_;

        PtrList<uniformDimensionedScalarField> level0EdgePtrs_;

        PtrList<refinementHistory> refHistoryPtrs_;


public:

    // Constructors

        //- Construct read. Assumes only first mesh in list
        //  having the files
        hexRef8DataList
        (
            PtrList<fvMesh>& myProcMeshes,
            const labelList& meshToProc,
            const labelListList& meshNoInProc
        );

        //- Disallow default bitwise copy construction
        hexRef8DataList(const hexRef8DataList&) = delete;


    //- Destructor
    ~hexRef8DataList();


    // Member Functions

        //- In-place distribute
        void distribute
        (
            PtrList<polyDistributionMap>& myProcDist
        );

        //- In-place reorder
        void reorder
        (
            PtrList<polyTopoChangeMap>& myProcMap
        );

        //- Write
        bool write() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const hexRef8DataList&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
