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
    (c) 2015-2019 OpenCFD Ltd.
    (c) 2022-2025 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::FaceZoneDataSetProvider

Description
    Provider for the face zone patches, which computes their vtk dataset

SourceFiles
    faceZoneDataSetProvider.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "itemDataSetProviders/mesh/meshDataSetProvider.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations

namespace Foam
{
class faceZone;
}

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                   Class FaceZoneDataSetProvider Declaration
\*---------------------------------------------------------------------------*/

class FaceZoneDataSetProvider : public MeshDataSetProvider
{
protected:

    // Looks up and interpolates a GeometricField from mesh_, then creates
    // and populates a field with the faceZoneMesh data, then converts that
    // field to vtk data and adds it as a piece to multiPieceDataSet_.
    // Templated type "fieldType" is one of scalar or vector.
    template<class fieldType>
    void addFaceZoneDataToConvertedMesh
    (
        const word& foamFieldName,
        const faceZone& myFaceZone,
        vtkDataSet* convertedMesh
    );

    void update(scalar currentTime) override;

public:
    FaceZoneDataSetProvider(const string& name, const string& region, const MeshAndFields& meshAndFields)
            : MeshDataSetProvider(name, region, meshAndFields) {};
    ~FaceZoneDataSetProvider() override = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
