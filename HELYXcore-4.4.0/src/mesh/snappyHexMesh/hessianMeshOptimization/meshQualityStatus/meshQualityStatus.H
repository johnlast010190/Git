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
    Foam::meshQualityStatus

Description
    Class that contains information for the current mesh quality used by mesh optimization
    algorithm

SourceFiles
    meshQualityStatus.C

\*---------------------------------------------------------------------------*/

#ifndef MESHQUALITYSTATUS_H_
#define MESHQUALITYSTATUS_H_

#include "fvMesh/fvMesh.H"
#include "cfdTools/general/include/fvCFD.H"
#include "motionSmoother/polyMeshGeometry/polyMeshGeometry.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "hessianMeshOptimization/meshMetric/meshMetric.H"

namespace Foam
{

class polyMeshGeometry;

/*---------------------------------------------------------------------------*\
                          Class meshQualityStatus Declaration
\*---------------------------------------------------------------------------*/
class meshQualityStatus
{
    const fvMesh& mesh_;

    const meshMetric& metric_;

    const dictionary& dict_;

    label nWrongFaces_;

    boolList isAWrongFace_;

    boolList activePoints_;

    // Private Member Function
    void findWrongFaces(labelHashSet& wrongFaces);

    void calculateActivePointSet();

    void incrementActivePointSet();

    tmp<vectorField> calculateCorrectCellCentres();

public  :

        meshQualityStatus
        (
            const fvMesh& mesh,
            const meshMetric& metric,
            const dictionary& dict
        );

        virtual ~meshQualityStatus();

   // Member Functions

        void checkMeshQuality();

        const boolList& getActivePoints() const;
        const boolList& getActiveFaces() const;
//        const boolList& getActiveCells() const;

        const boolList& getWrongFaces() const;
        const label& numberOfWrongFaces() const;

};


} /* namespace Foam */

#endif /* MESHQUALITYSTATUS_H_ */
