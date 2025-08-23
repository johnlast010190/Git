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
    (c) 2025 OpenCFD Ltd.

Class
    Foam::functionObjects::runTimeVis::vtuAdaptor

Description
    A low-level backend adaptor for converting OpenFOAM volume meshes/fields
    to VTK internal representation.
    The output is a single piece vtkUnstructuredGrid dataset.
    Unlike foamVtkVtuAdaptor, the points are wrapped instead of copied

SourceFiles
    rtppVtkVtuAdaptor.C

\*---------------------------------------------------------------------------*/

#ifndef rtpp_vtuAdaptor_H
#define rtpp_vtuAdaptor_H

#include "fvMesh/fvMesh.H"
#include "vtk/adaptor/foamVtkTools.H"
#include "vtk/part/foamVtkMeshMaps.H"
#include "vtk/adaptor/foamVtkVtuAdaptor.H"

#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkUnstructuredGrid.h"

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

class vtkCellArray;
class vtkDataSet;
class vtkFloatArray;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                       Class vtk::vtuAdaptor Declaration
\*---------------------------------------------------------------------------*/

//- Bookkeeping for vtkUnstructuredGrid
class rtppVtuAdaptor
{
public:

    explicit rtppVtuAdaptor(Foam::vtk::vtuAdaptor& adaptor) : adaptor(adaptor)
    {};
    // Member Functions

    // Mesh Conversion

        //- The vtk points for the mesh (and decomposition)
        vtkSmartPointer<vtkPoints> points
        (
            const fvMesh& mesh
        ) const;

        //- Internal mesh as vtkUnstructuredGrid.
        //  Calling this automatically populates the sizing and vtk::Caching
        vtkSmartPointer<vtkUnstructuredGrid> internal
        (
            const fvMesh& mesh,
            bool decompPoly = false
        ) const;

protected:
    Foam::vtk::vtuAdaptor& adaptor;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
