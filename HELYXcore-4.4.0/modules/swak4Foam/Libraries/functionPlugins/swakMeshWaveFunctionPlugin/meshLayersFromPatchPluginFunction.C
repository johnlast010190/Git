/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) ICE Stroemungsfoschungs GmbH
    (c) 1991-2008 OpenCFD Ltd.

Contributors/Copyright:
    2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "meshLayersFromPatchPluginFunction.H"
#include "FieldValueExpressionDriver.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

namespace Foam {

defineTypeNameAndDebug(meshLayersFromPatchPluginFunction,1);
addNamedToRunTimeSelectionTable(FieldValuePluginFunction, meshLayersFromPatchPluginFunction , name, meshLayersFromPatch);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshLayersFromPatchPluginFunction::meshLayersFromPatchPluginFunction(
    const FieldValueExpressionDriver &parentDriver,
    const word &name
):
    meshLayersGeneralPluginFunction(
        parentDriver,
        name,
        string("patchName primitive word")
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void meshLayersFromPatchPluginFunction::setArgument(
    label index,
    const word &patchName
) {
    assert(index==0);

    patchName_=patchName;
}

void meshLayersFromPatchPluginFunction::initFacesAndCells()
{
    label patchI=mesh().boundaryMesh().findPatchID(patchName_);
    if (patchI<0) {
        FatalErrorIn("meshLayersFromPatchPluginFunction::initFacesAndCells()")
            << "Patch name " << patchName_ << " not in valid names"
                << mesh().boundaryMesh().names()
                << endl
                << exit(FatalError);
    }

    startFaces_=labelList(mesh().boundaryMesh()[patchI].size());
    for (label i=0;i<mesh().boundaryMesh()[patchI].size();i++) {
        startFaces_[i]=mesh().boundaryMesh()[patchI].start()+i;
    }
    startValues_=List<MeshLayersDistFromPatch>(
        mesh().boundaryMesh()[patchI].size(),
        MeshLayersDistFromPatch(1)
    );
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace

// ************************************************************************* //
