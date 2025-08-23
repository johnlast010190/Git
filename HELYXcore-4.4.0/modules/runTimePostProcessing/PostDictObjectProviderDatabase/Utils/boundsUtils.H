/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.0.1
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

namespace
    Foam::functionObjects::runTimeVis::boundsUtils

Description
    Class for calculating parameters to do with dataset bounds

SourceFiles
    boundsUtils.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

#include "primitives/Scalar/scalar/scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class vtkDataSet;
class vtkRenderer;

namespace Foam
{

class boundBox;

/*---------------------------------------------------------------------------*\
                           boundsUtils Declaration
\*---------------------------------------------------------------------------*/
namespace functionObjects::runTimeVis::boundsUtils
{
    // Calculates the bounds for a polyData considering all processors
    void getFullPolyBounds(vtkDataSet* poly, scalar result[6]);

    boundBox computeRenderBoundingBox(vtkRenderer* renderer);
    boundBox computeAllProcsRenderBoundingBox(vtkRenderer* renderer);

    void computeRenderBoundingBox(vtkRenderer* renderer, scalar dBounds[6]);
    void computeAllProcsRenderBoundingBox(vtkRenderer* renderer, scalar dBounds[6]);

    boundBox boundBoxFromDoubleArray(double dBounds[6]);
    void doubleArrayFromBoundBox(double dBounds[6], const boundBox& bounds);
    void scalarArrayFromBoundBox(scalar sBounds[6], const boundBox& bounds);
}

} // End namespace
