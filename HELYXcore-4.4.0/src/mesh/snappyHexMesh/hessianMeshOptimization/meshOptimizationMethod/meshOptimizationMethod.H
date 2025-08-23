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
    (c) 2014 Engys Ltd.

Class
    Foam::meshOptimizationMethod

Description
    Abstract class for choosing the mesh optimization method

SourceFiles
    meshOptimizationMethod.C

\*---------------------------------------------------------------------------*/

#ifndef meshOptimizationMethod_H
#define meshOptimizationMethod_H

#include "fvMesh/fvMesh.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "hessianMeshOptimization/meshOptimizationMethod/activeSet/activeSet.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "hessianMeshOptimization/meshMetric/meshMetric.H"

namespace Foam
{
/*---------------------------------------------------------------------------*\
                   Class meshOptimizationMethod Declaration
\*---------------------------------------------------------------------------*/
class meshOptimizationMethod
{

    // Private data

        //- Mesh
        const fvMesh& mesh_;


protected:

        activeSet activeSet_;

        //Given a mesh boolList shrink it to keep only the true labels
        labelList compactList(const boolList& bL);

        //Increments the activePoints by the cellPoints of each active point
        void incrementActivePointsSet(boolList& activePoints);


public:

    //- Runtime type information
    TypeName("meshOptimizationMethod");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        meshOptimizationMethod,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const meshMetric& metric
        ),
        (dict,mesh, metric)
    );


    // Constructors

        //- Construct from components
        meshOptimizationMethod
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const meshMetric& metric
        );

        //- Disallow default bitwise copy construction
        meshOptimizationMethod(const meshOptimizationMethod&) = delete;


    //- Destructor
    virtual ~meshOptimizationMethod();


    // Selectors

        //- Select null constructed
        static autoPtr<meshOptimizationMethod> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const meshMetric& metric
        );

        static autoPtr<meshOptimizationMethod> NewFull
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const meshMetric& metric
        );

    // Member Functions

        virtual bool updateActiveSet() = 0;

        virtual bool converge() = 0;

        virtual void merge(const labelList& pointSet) = 0;

        const activeSet& getActiveSet() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshOptimizationMethod&) = delete;
};

} /* namespace Foam */

#endif
