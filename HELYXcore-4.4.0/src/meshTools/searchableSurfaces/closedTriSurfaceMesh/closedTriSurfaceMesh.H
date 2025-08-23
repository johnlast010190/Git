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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::closedTriSurfaceMesh

Description
    A triSurfaceMesh where it is forced to check volumeTypes, used for surfaces
    that are topologically non-manifold (small holes or multiple parts) but are
    geometrically essentially closed

SourceFiles
    closedTriSurfaceMesh.C

\*---------------------------------------------------------------------------*/

#ifndef closedTriSurfaceMesh_H
#define closedTriSurfaceMesh_H

#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class closedTriSurfaceMesh Declaration
\*---------------------------------------------------------------------------*/

class closedTriSurfaceMesh
:
    public triSurfaceMesh
{
public:

    //- Runtime type information
    TypeName("closedTriSurfaceMesh");


    // Constructors

        //- Construct from triSurface
        closedTriSurfaceMesh(const IOobject&, const triSurface&);

        //- Construct read.
        closedTriSurfaceMesh(const IOobject& io);

        //- Construct from IO and dictionary (used by searchableSurface).
        //  Dictionary may contain a 'scale' entry (eg, 0.001: mm -> m)
        closedTriSurfaceMesh
        (
            const IOobject& io,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        closedTriSurfaceMesh(const closedTriSurfaceMesh&) = delete;


    //- Destructor
    virtual ~closedTriSurfaceMesh();


    // Member Functions

        //- Whether supports volume type, forcing to true to force getVolumeType
        //  queries for this type
        virtual bool hasVolumeType() const
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const closedTriSurfaceMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
