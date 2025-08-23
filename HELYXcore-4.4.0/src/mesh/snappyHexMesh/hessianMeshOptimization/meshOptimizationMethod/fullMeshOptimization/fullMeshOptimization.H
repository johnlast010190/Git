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
    (c) 2014-2024 Engys Ltd.

Class
    Foam::fullmeshOptimization

Description
    Optimize every point in the fvMesh

SourceFiles
    meshOptimizationMethod.C

\*---------------------------------------------------------------------------*/
#ifndef fullMeshOptimization_H
#define fullMeshOptimization_H

#include "fvMesh/fvMesh.H"
#include "hessianMeshOptimization/meshOptimizationMethod/meshOptimizationMethod.H"

namespace Foam
{
/*---------------------------------------------------------------------------*\
                   Class fullMeshOptimization Declaration
\*---------------------------------------------------------------------------*/
class fullMeshOptimization
:
    public meshOptimizationMethod
{
    // Private Data

        const fvMesh& mesh_;

        void calculateActiveSet();


public:

    //- Runtime type information
    TypeName("fullMeshOptimization");


    // Constructors

        //- Construct from components
        fullMeshOptimization
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const meshMetric& metric
        );

        //- Disallow default bitwise copy construction
        fullMeshOptimization(const fullMeshOptimization&) = delete;


    //- Destructor
    virtual ~fullMeshOptimization();


    // Public Member Functions

        bool updateActiveSet();

        bool converge();

        void merge(const labelList& pointSet);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fullMeshOptimization&) = delete;
};

} /* namespace Foam */

#endif /*fullMeshOptimization_H */
