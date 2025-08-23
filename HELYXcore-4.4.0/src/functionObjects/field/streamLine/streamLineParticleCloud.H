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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::streamLineParticleCloud

Description
    A Cloud of streamLine particles

SourceFiles
    streamLineCloud.C

\*---------------------------------------------------------------------------*/

#ifndef streamLineCloud_H
#define streamLineCloud_H

#include "Cloud/Cloud.H"
#include "streamLine/streamLineParticle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class streamLineCloud Declaration
\*---------------------------------------------------------------------------*/

class streamLineParticleCloud
:
    public Cloud<streamLineParticle>
{
public:

    //- Type of parcel the cloud was instantiated for
    typedef streamLineParticle parcelType;


    // Constructors

        //- Construct given mesh
        streamLineParticleCloud
        (
            const polyMesh&,
            const word& cloudName = "defaultCloud",
            bool readFields = true
        );

        //- Construct from mesh, cloud name, and a list of particles
        streamLineParticleCloud
        (
            const polyMesh& mesh,
            const word& cloudName,
            const IDLList<streamLineParticle>& particles
        );

        //- Disallow default bitwise copy construction
        streamLineParticleCloud(const streamLineParticleCloud&) = delete;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const streamLineParticleCloud&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
