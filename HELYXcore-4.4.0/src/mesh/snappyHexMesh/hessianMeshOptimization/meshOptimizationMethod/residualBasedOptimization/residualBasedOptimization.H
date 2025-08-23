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
    Foam::residualBasedOptimization

Description
    Optimizing the points that are actually moving during the optimization

SourceFiles
    qualityBasedOptimization.C

\*---------------------------------------------------------------------------*/

#ifndef residualBasedOptimization_H
#define qualityBasedOptimization_H

#include "fvMesh/fvMesh.H"
#include "hessianMeshOptimization/meshOptimizationMethod/meshOptimizationMethod.H"
#include "hessianMeshOptimization/meshQualityStatus/meshQualityStatus.H"

namespace Foam
{
class residualBasedOptimization
:
    public meshOptimizationMethod
{
    // Private Data

    const fvMesh& mesh_;

    const pointField& points_;

    pointField oldPoints_;

    scalarField maxDisp_;

    labelField counter_;

    label pL_;

    void initialActiveSet();

    void calculateActiveSet();


public:

    //- Runtime type information
    TypeName("residualBasedOptimization");

    // Constructors

        //- Construct from a dictionary, fvMesh and pointField
        residualBasedOptimization
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const pointField& points
        );

        //- Disallow default bitwise copy construction
        residualBasedOptimization(const residualBasedOptimization&) = delete;


    //- Destructor
    virtual ~residualBasedOptimization();


    // Public Member Functions

        bool updateActiveSet();

        bool converge();

        void merge(const labelList& pointSet);


    // Member Operators

    //- Disallow default bitwise assignment
    void operator=(const residualBasedOptimization&) = delete;
};

} /* namespace Foam */

#endif /* residualBasedOptimization_H */
