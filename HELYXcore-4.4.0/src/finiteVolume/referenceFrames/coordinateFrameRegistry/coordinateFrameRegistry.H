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
    (c) 2020-2024 Engys Ltd.

Class
    Foam::coordinateFrameRegistry

Description

\*---------------------------------------------------------------------------*/

#ifndef coordinateFrameRegistry_H
#define coordinateFrameRegistry_H

#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class coordinateFrameRegistry Declaration
\*---------------------------------------------------------------------------*/

class coordinateFrameRegistry
{
protected:

    // Protected Member Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- List of regIOobjects name that use this frame
        mutable DynamicList<word> registeredNames_;

        //- List of patches attached to this frame
        mutable DynamicList<label> registeredPatches_;

        //- Flag to update the faces
        mutable bool faceSetNeedUpdate_;


public:


    // Constructors

        //- Construct from mesh
        coordinateFrameRegistry(const fvMesh& mesh);


    //- Destructor
    virtual ~coordinateFrameRegistry() = default;


    // Member Functions

        //- Registered names
        const DynamicList<word>& registeredNames() const
        {
            return registeredNames_;
        }

        //- Attach object to the frame
        void attachObject(const word& name) const;

        //- Attach patch to the frame
        bool attachPatch(const label patchi) const;

        //- Check if the patch is attached to the MRF/GRF
        bool isAttachToMRF(const label patchi) const;

        //- Update lists of faces
        void updateFaces() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif
