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

#include "floodFillFromFacesPluginFunction.H"
#include "FieldValueExpressionDriver.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchFields.H"

namespace Foam {

defineTypeNameAndDebug(floodFillFromFacesPluginFunction,1);
addNamedToRunTimeSelectionTable(FieldValuePluginFunction, floodFillFromFacesPluginFunction , name, floodFillFromFaces);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

floodFillFromFacesPluginFunction::floodFillFromFacesPluginFunction(
    const FieldValueExpressionDriver &parentDriver,
    const word &name
):
    floodFillGeneralPluginFunction(
        parentDriver,
        name,
        string("blockedFaces internalField surfaceLogicalField")
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void floodFillFromFacesPluginFunction::setArgument(
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

void floodFillFromFacesPluginFunction::initFacesAndCells()
{
    labelHashSet facesBlocked;
    const cellList &cells=mesh().cells();

    forAll(blocked_(),faceI) {
        if (CommonValueExpressionDriver::toBool(blocked_()[faceI])) {
            facesBlocked.insert(faceI);
            faceValues_[faceI]=FloodFillData(0,0,true);
        }
    }
    forAll(blocked_().boundaryField(),patchI) {
        const fvsPatchField<scalar> &p=blocked_().boundaryField()[patchI];
        forAll(p,i) {
            if (CommonValueExpressionDriver::toBool(p[i])) {
                label faceI=i+p.patch().patch().start();
                facesBlocked.insert(faceI);
                faceValues_[faceI]=FloodFillData(0,0,true);
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
            cellValues_[cellI]=FloodFillData(0,0,true);
        }
    }
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace

// ************************************************************************* //
