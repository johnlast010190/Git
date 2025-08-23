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

\*---------------------------------------------------------------------------*/

#include "boundsUtils.H"

#include "primitives/Vector/Vector.H"
#include "meshes/boundBox/boundBox.H"
#include "Utils/ParallelUtils.H"

#include "vtkRenderer.h"
#include "vtkActorCollection.h"

namespace Foam::functionObjects::runTimeVis
{

void boundsUtils::getFullPolyBounds(vtkDataSet* poly, scalar result[6]) {
    if (ParallelUtils::isRunningInParallel()) {
        // Can't use Foam::boundBox here, because that would truncate our
        // doubles to floats when in SP mode.  Perhaps we should template
        // boundBox?
        double bounds[6];
        if (poly->GetNumberOfPoints() > 0)
        {
            poly->GetBounds(bounds);
        }
        else
        {
            for (int i = 0; i < 6; i += 2)
            {
                bounds[i] = VTK_DOUBLE_MAX;
                bounds[i+1] = VTK_DOUBLE_MIN;
            }
        }

        Vector<scalar> minPoint = { scalar(bounds[0]), scalar(bounds[2]), scalar(bounds[4]) };
        Vector<scalar> maxPoint = { scalar(bounds[1]), scalar(bounds[3]), scalar(bounds[5]) };
        reduce(minPoint, minOp<Vector<scalar>>());
        reduce(maxPoint, maxOp<Vector<scalar>>());
        // Not superbly elegant...
        result[0] = minPoint[0];
        result[2] = minPoint[1];
        result[4] = minPoint[2];
        result[1] = maxPoint[0];
        result[3] = maxPoint[1];
        result[5] = maxPoint[2];
    } else {
        double bounds[6];
        poly->GetBounds(bounds);
        for (int i = 0; i < 6; i++) result[i] = static_cast<scalar>(bounds[i]);
    }
}

boundBox boundsUtils::computeRenderBoundingBox(vtkRenderer* renderer)
{
    boundBox bounds;
    double dBounds[6];
    vtkActorCollection *actors = renderer->GetActors();
    actors->InitTraversal();
    // get the bounds for each actor
    for (vtkActor *actor = actors->GetNextActor(); actor != nullptr; actor = actors->GetNextActor())
    {
        actor->GetBounds(dBounds);
        bounds.add(boundBoxFromDoubleArray(dBounds));
    }
    return bounds;
}

void boundsUtils::computeRenderBoundingBox(vtkRenderer* renderer, scalar dBounds[6])
{
    scalarArrayFromBoundBox(dBounds, computeRenderBoundingBox(renderer));
}

boundBox boundsUtils::computeAllProcsRenderBoundingBox(vtkRenderer* renderer)
{
    boundBox bounds = computeRenderBoundingBox(renderer);
    bounds.reduce();
    return bounds;
}

void boundsUtils::computeAllProcsRenderBoundingBox(vtkRenderer* renderer, scalar dBounds[6])
{
    scalarArrayFromBoundBox(dBounds, computeAllProcsRenderBoundingBox(renderer));
}



boundBox boundsUtils::boundBoxFromDoubleArray(double dBounds[6])
{
    if(dBounds[0] > dBounds[1] || dBounds[2] > dBounds[3] || dBounds[4] > dBounds[5])
    {
        return {};
    }
    return {point(static_cast<scalar>(dBounds[0]),
                  static_cast<scalar>(dBounds[2]),
                  static_cast<scalar>(dBounds[4])),
            point(static_cast<scalar>(dBounds[1]),
                  static_cast<scalar>(dBounds[3]),
                  static_cast<scalar>(dBounds[5]))
    };
}

void boundsUtils::doubleArrayFromBoundBox(double dBounds[6], const boundBox& bounds)
{
    dBounds[0] = bounds.min().x();
    dBounds[1] = bounds.max().x();
    dBounds[2] = bounds.min().y();
    dBounds[3] = bounds.max().y();
    dBounds[4] = bounds.min().z();
    dBounds[5] = bounds.max().z();
}

void boundsUtils::scalarArrayFromBoundBox(scalar sBounds[6], const boundBox& bounds)
{
    sBounds[0] = bounds.min().x();
    sBounds[1] = bounds.max().x();
    sBounds[2] = bounds.min().y();
    sBounds[3] = bounds.max().y();
    sBounds[4] = bounds.min().z();
    sBounds[5] = bounds.max().z();
}

} // End namespace
