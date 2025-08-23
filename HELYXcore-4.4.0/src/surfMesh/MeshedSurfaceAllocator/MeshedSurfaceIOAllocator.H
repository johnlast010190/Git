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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::MeshedSurfaceIOAllocator

Description
    A helper class for storing points, faces and zones with IO capabilities.

SourceFiles
    MeshedSurfaceIOAllocator.C

\*---------------------------------------------------------------------------*/

#ifndef MeshedSurfaceIOAllocator_H
#define MeshedSurfaceIOAllocator_H

#include "meshes/primitiveShapes/point/pointIOField.H"
#include "meshes/meshShapes/face/faceIOList.H"
#include "surfZone/surfZone/surfZoneIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class MeshedSurfaceIOAllocator Declaration
\*---------------------------------------------------------------------------*/

class MeshedSurfaceIOAllocator
{
    // Private data

        //- Points
        pointIOField points_;

        //- Faces
        faceCompactIOList faces_;

        //- Surface zones
        surfZoneIOList zones_;


public:

    // Constructors

        //- Read construct from IOobjects
        MeshedSurfaceIOAllocator
        (
            const IOobject& ioPoints,
            const IOobject& ioFaces,
            const IOobject& ioZones
        );

        //- Construct from IOobjects, copying components
        MeshedSurfaceIOAllocator
        (
            const IOobject& ioPoints,
            const pointField& points,
            const IOobject& ioFaces,
            const faceList& faces,
            const IOobject& ioZones,
            const surfZoneList& zones
        );

        //- Construct from IOobjects, possibly transferring components
        MeshedSurfaceIOAllocator
        (
            const IOobject& ioPoints,
            pointField&& points,
            const IOobject& ioFaces,
            faceList&& faces,
            const IOobject& ioZones,
            surfZoneList&& zones
        );

        //- Disallow default bitwise copy construction
        MeshedSurfaceIOAllocator(const MeshedSurfaceIOAllocator&) = delete;


    //- Destructor
    virtual ~MeshedSurfaceIOAllocator();


    // Member Functions

        // Database

            //- Set the instance for mesh files
            void setInstance(const fileName&);

            //- Adjust the write option for all components
            void setWriteOption(IOobject::writeOption);


        // Access

            //- Non-const access to the points
            pointIOField& storedIOPoints()
            {
                return points_;
            }

            //- Non-const access to the faces
            faceCompactIOList& storedIOFaces()
            {
                return faces_;
            }

            //- Non-const access to the zones
            surfZoneIOList& storedIOZones()
            {
                return zones_;
            }

            //- Const access to the points
            const pointIOField& storedIOPoints() const
            {
                return points_;
            }

            //- Const access to the faces
            const faceCompactIOList& storedIOFaces() const
            {
                return faces_;
            }

            //- Const access to the zones
            const surfZoneIOList& storedIOZones() const
            {
                return zones_;
            }


        //  Storage management

            //- Clear primitive data (points, faces and zones)
            void clear();

            //- Reset primitive data (faces and zones)
            //  Note, optimized to avoid overwriting data
            void resetFaces
            (
                faceList&& faces,
                surfZoneList&& zones
            );

            //- Reset primitive data (points, faces and zones)
            //  Note, optimized to avoid overwriting data
            void reset
            (
                pointField&& points,
                faceList&& faces,
                surfZoneList&& zones
            );

            //- Reset primitive data (points, faces and zones)
            //  Note, optimized to avoid overwriting data
            void reset
            (
                List<point>&& points,
                faceList&& faces,
                surfZoneList&& zones
            );


        // Writing

            //- Write using given format, version and compression
            virtual bool writeObject
            (
                IOstream::streamFormat,
                IOstream::versionNumber,
                IOstream::compressionType,
                const bool write
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const MeshedSurfaceIOAllocator&) = delete;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
