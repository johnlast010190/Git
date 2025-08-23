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
    (c) 2018-2024 Engys Ltd.
    (c) 2016-2017 Wikki Ltd

Class
    Foam::faGlobalMeshData

Description
    Various mesh related information for a parallel run.

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#ifndef faGlobalMeshData_H
#define faGlobalMeshData_H

#include "faMesh/faGlobalMeshData/faProcessorTopology.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class faMesh;

/*---------------------------------------------------------------------------*\
                        Class faGlobalMeshData Declaration
\*---------------------------------------------------------------------------*/

class faGlobalMeshData
:
    public faProcessorTopology
{
    // Private data

        //- Reference to mesh
        const faMesh& mesh_;

        // Globally shared point addressing

            //- Total number of global points
            label nGlobalPoints_;

            //- Indices of local points that are globally shared
            labelList sharedPointLabels_;

            //- Indices of globally shared points in the master list
            //  This list contains all the shared points in the mesh
            labelList sharedPointAddr_;


public:

    //- Runtime type information
    ClassName("faGlobalMeshData");


    // Constructors

        //- Construct from mesh
        faGlobalMeshData(const faMesh& mesh);

        //- Disallow default bitwise copy construction
        faGlobalMeshData(const faGlobalMeshData&) = delete;


    //- Destructor
    ~faGlobalMeshData();


    // Member Functions

        // Access

            //- Return mesh reference
            const faMesh& mesh() const;

        // Globally shared point addressing

            //- Return number of globally shared points
            label nGlobalPoints() const
            {
                return nGlobalPoints_;
            }

            //- Return indices of local points that are globally shared
            const labelList& sharedPointLabels() const
            {
                return sharedPointLabels_;
            }

            //- Return addressing into the complete globally shared points
            //  list
            const labelList& sharedPointAddr() const
            {
                return sharedPointAddr_;
            }

            //- Change global mesh data given a topological change.
            void topoChange();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const faGlobalMeshData&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
