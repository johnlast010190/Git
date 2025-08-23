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
    (c) 2014-2016 OpenFOAM Foundation

Class
    Foam::decompositionModel

Description
    MeshObject wrapper of decompositionMethod

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef decompositionModel_H
#define decompositionModel_H

#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "meshes/MeshObject/MeshObject.H"
#include "decompositionMethod/decompositionMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyTopoChangeMap;
class polyMesh;

/*---------------------------------------------------------------------------*\
                      Class decompositionModel Declaration
\*---------------------------------------------------------------------------*/

class decompositionModel
:
    public MeshObject
    <
        polyMesh,
        UpdateableMeshObject,
        decompositionModel
    >,
    public IOdictionary
{

    // Private data

        mutable autoPtr<decompositionMethod> decomposerPtr_;


public:

    // Declare name of the class and its debug switch
    ClassName("decompositionModel");


    // Selectors

        //- Read (optionallly from absolute path) & register on mesh
        static const decompositionModel& New
        (
            const polyMesh& mesh,
            const fileName& decompDictFile = ""
        );

        //- Read (optionallly from supplied dictionary) & register on mesh
        static const decompositionModel& New
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const fileName& decompDictFile = ""
        );


    // Constructors

        //- Construct from typeName or optional path to controlDictionary
        decompositionModel(const polyMesh&, const fileName& = "");


        //- Construct from typeName or optional path to controlDictionary
        decompositionModel
        (
            const polyMesh&,
            const dictionary& dict,
            const fileName& = ""
        );


    // Member functions

        using IOdictionary::found;

        decompositionMethod& decomposer() const
        {
            if (!decomposerPtr_.valid())
            {
                decomposerPtr_ = decompositionMethod::New(*this);
            }
            return decomposerPtr_();
        }

        //- Helper: return IOobject with optionally absolute path provided
        static IOobject selectIO(const IOobject&, const fileName&);


        // UpdateableMeshObject

            virtual bool movePoints()
            {
                return false;
            }

            virtual void topoChange(const polyTopoChangeMap&)
            {}

            virtual void mapMesh(const polyMeshMap&)
            {}

            virtual void distribute(const polyDistributionMap&)
            {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
