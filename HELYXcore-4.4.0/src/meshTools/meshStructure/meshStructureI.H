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
    (c) 2013 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "meshStructure/meshStructure.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::meshStructure::structured() const
{
    return structured_;
}


const Foam::labelList& Foam::meshStructure::cellToPatchFaceAddressing() const
{
    return cellToPatchFaceAddressing_;
}


Foam::labelList& Foam::meshStructure::cellToPatchFaceAddressing()
{
    return cellToPatchFaceAddressing_;
}


const Foam::labelList& Foam::meshStructure::cellLayer() const
{
    return cellLayer_;
}


Foam::labelList& Foam::meshStructure::cellLayer()
{
    return cellLayer_;
}


const Foam::labelList& Foam::meshStructure::faceToPatchFaceAddressing() const
{
    return faceToPatchFaceAddressing_;
}


Foam::labelList& Foam::meshStructure::faceToPatchFaceAddressing()
{
    return faceToPatchFaceAddressing_;
}


const Foam::labelList& Foam::meshStructure::faceToPatchEdgeAddressing() const
{
    return faceToPatchEdgeAddressing_;
}


Foam::labelList& Foam::meshStructure::faceToPatchEdgeAddressing()
{
    return faceToPatchEdgeAddressing_;
}


const Foam::labelList& Foam::meshStructure::faceLayer() const
{
    return faceLayer_;
}


Foam::labelList& Foam::meshStructure::faceLayer()
{
    return faceLayer_;
}


const Foam::labelList& Foam::meshStructure::pointToPatchPointAddressing() const
{
    return pointToPatchPointAddressing_;
}


Foam::labelList& Foam::meshStructure::pointToPatchPointAddressing()
{
    return pointToPatchPointAddressing_;
}


const Foam::labelList& Foam::meshStructure::pointLayer() const
{
    return pointLayer_;
}


Foam::labelList& Foam::meshStructure::pointLayer()
{
    return pointLayer_;
}


// ************************************************************************* //
