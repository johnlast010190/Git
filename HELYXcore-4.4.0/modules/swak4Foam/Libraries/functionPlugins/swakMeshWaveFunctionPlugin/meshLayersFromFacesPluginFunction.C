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

#include "meshLayersFromFacesPluginFunction.H"
#include "FieldValueExpressionDriver.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchFields.H"

namespace Foam {

defineTypeNameAndDebug(meshLayersFromFacesPluginFunction,1);
addNamedToRunTimeSelectionTable(FieldValuePluginFunction, meshLayersFromFacesPluginFunction , name, meshLayersFromFaces);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshLayersFromFacesPluginFunction::meshLayersFromFacesPluginFunction(
    const FieldValueExpressionDriver &parentDriver,
    const word &name
):
    meshLayersGeneralPluginFunction(
        parentDriver,
        name,
        string("blockedFaces internalField surfaceLogicalField")
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void meshLayersFromFacesPluginFunction::setArgument(
    label index,
    const string &content,
    const CommonValueExpressionDriver &driver
) {
    assert(index==0);

    blocked_.set(
        new surfaceScalarField(
            dynamic_cast<const FieldValueExpressionDriver &>(
                driver
            ).getResult<surfaceScalarField>()
        )
    );
}

void meshLayersFromFacesPluginFunction::initFacesAndCells()
{
    labelHashSet facesBlocked;
    const cellList &cells=mesh().cells();

    forAll(blocked_(),faceI) {
        if (CommonValueExpressionDriver::toBool(blocked_()[faceI])) {
            facesBlocked.insert(faceI);
        }
    }
    forAll(blocked_().boundaryField(),patchI) {
        const fvsPatchField<scalar> &p=blocked_().boundaryField()[patchI];
        forAll(p,i) {
            if (CommonValueExpressionDriver::toBool(p[i])) {
                facesBlocked.insert(i+p.patch().patch().start());
            }
        }
    }
    forAll(cells,cellI) {
        bool allBlocked=true;
        forAll(cells[cellI],i) {
            if (!facesBlocked.found(cells[cellI][i])) {
                allBlocked=false;
                break;
            }
        }
        if (allBlocked) {
            cellValues_[cellI]=MeshLayersDistFromPatch(0,true);
        }
    }

    startFaces_=labelList(facesBlocked.toc());
    forAll(startFaces_,i) {
        faceValues_[startFaces_[i]]=MeshLayersDistFromPatch(0,true);
    }
    startValues_=List<MeshLayersDistFromPatch>(
        startFaces_.size(),
        MeshLayersDistFromPatch(1,false)
    );
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace

// ************************************************************************* //
