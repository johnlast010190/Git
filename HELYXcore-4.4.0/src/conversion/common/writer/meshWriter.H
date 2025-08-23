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
    (c) 2016 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::meshWriter

Description
    Write OpenFOAM meshes and/or results to another CFD format.

\par Files

    "constant/boundaryRegion" is an IOMap<dictionary> that contains
    the boundary type and names. eg,
    \verbatim
        (
            0
            {
                BoundaryType    wall;
                Label           Default_Boundary_Region;
            }

            1
            {
                BoundaryType    inlet;
                Label           inlet_1;
            }

            ...

            4
            {
                BoundaryType    pressure;
                Label           outlet;
            }
        )
    \endverbatim


SourceFiles
    meshWriter.C

\*---------------------------------------------------------------------------*/

#ifndef meshWriter_H
#define meshWriter_H

#include "meshes/polyMesh/polyMesh.H"
#include "common/tables/boundaryRegion.H"
#include "common/tables/cellTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class meshWriter Declaration
\*---------------------------------------------------------------------------*/

class meshWriter
{
protected:

    // Protected data

        //- Mesh reference
        const polyMesh& mesh_;

        //- Scaling factor for points (eg, [m] -> [mm])
        scalar scaleFactor_;

        //- boundaryRegion persistent data saved as a dictionary
        boundaryRegion boundaryRegion_;

        //- cellTable persistent data saved as a dictionary
        cellTable cellTable_;

        //- cellTable IDs for each cell
        labelList cellTableId_;

        //- Pointers to cell shape models
        static const cellModel* unknownModel;
        static const cellModel* tetModel;
        static const cellModel* pyrModel;
        static const cellModel* prismModel;
        static const cellModel* hexModel;


public:

    // Static data members

        //- Specify a default mesh name
        static string defaultMeshName;


    // Constructors

        //- Create a writer object with given output scaling
        meshWriter
        (
            const polyMesh&,
            const scalar scaling = 1.0
        );

        //- Disallow default bitwise copy construction
        meshWriter(const meshWriter&) = delete;


    //- Destructor
    virtual ~meshWriter();


    // Member Functions

        //- Write volume mesh. Subclass must supply this method
        virtual bool write
        (
            const fileName& timeName = fileName::null
        ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshWriter&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
