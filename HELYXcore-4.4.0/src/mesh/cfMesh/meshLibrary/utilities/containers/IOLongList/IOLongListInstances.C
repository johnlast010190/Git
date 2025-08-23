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
    (c) Creative Fields, Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Description
    Declaration of IOLongList ClassNames for IOLists that do not have .C files.

\*---------------------------------------------------------------------------*/

#include "utilities/containers/IOLongList/IOLongListInstances.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineCompoundTypeName(IOLongList<label>, labelIOListPMG);
    defineCompoundTypeName(IOLongList<point>, pointIOFieldPMG);
    //defineCompoundTypeName(IOLongList<face>, faceIOListPMG);
    //defineCompoundTypeName(IOLongList<cell>, cellIOListPMG);
    //addCompoundToRunTimeSelectionTable(IOLongList<label>, labelIOLongList);

    defineTemplateTypeNameAndDebugWithName(labelIOListPMG, "labelList", 0);
    defineTemplateTypeNameAndDebugWithName(pointIOFieldPMG, "vectorField", 0);
    //defineTemplateTypeNameAndDebugWithName(faceIOListPMG, "faceList", 0);
    //defineTemplateTypeNameAndDebugWithName(cellIOListPMG, "cellList", 0);
}

// ************************************************************************* //
