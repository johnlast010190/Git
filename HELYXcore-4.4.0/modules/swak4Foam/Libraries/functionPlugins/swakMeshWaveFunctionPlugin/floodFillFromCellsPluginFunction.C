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

#include "floodFillFromCellsPluginFunction.H"
#include "FieldValueExpressionDriver.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

namespace Foam {

defineTypeNameAndDebug(floodFillFromCellsPluginFunction,1);
addNamedToRunTimeSelectionTable(FieldValuePluginFunction, floodFillFromCellsPluginFunction , name, floodFillFromCells);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

floodFillFromCellsPluginFunction::floodFillFromCellsPluginFunction(
    const FieldValueExpressionDriver &parentDriver,
    const word &name
):
    floodFillGeneralPluginFunction(
        parentDriver,
        name,
        string("blockedCells internalField volLogicalField")
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void floodFillFromCellsPluginFunction::setArgument(
    label index,
    const string &content,
    const CommonValueExpressionDriver &driver
) {
    assert(index==0);

    blocked_.set(
        new volScalarField(
            dynamic_cast<const FieldValueExpressionDriver &>(
                driver
            ).getResult<volScalarField>()
        )
    );
}

void floodFillFromCellsPluginFunction::initFacesAndCells()
{
    labelHashSet facesBlocked;
    const cellList &cells=mesh().cells();

    forAll(blocked_(),cellI) {
        if (CommonValueExpressionDriver::toBool(blocked_()[cellI])) {
            cellValues_[cellI]=FloodFillData(0,0,true);
            const cell &c=cells[cellI];
            forAll(c,i) {
                facesBlocked.insert(c[i]);
            }
        }
    }

    forAll(startFaces_,i) {
        faceValues_[startFaces_[i]]=FloodFillData(0,0,true);
    }
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace

// ************************************************************************* //
