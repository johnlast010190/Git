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
    2012-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "mqCellShapePluginFunction.H"
#include "FieldValueExpressionDriver.H"

#include "meshes/meshShapes/cellMatcher/hexMatcher.H"
#include "meshes/meshShapes/cellMatcher/wedgeMatcher.H"
#include "meshes/meshShapes/cellMatcher/prismMatcher.H"
#include "meshes/meshShapes/cellMatcher/pyrMatcher.H"
#include "meshes/meshShapes/cellMatcher/tetWedgeMatcher.H"
#include "meshes/meshShapes/cellMatcher/tetMatcher.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

namespace Foam {

defineTypeNameAndDebug(mqCellShapePluginFunction,1);
addNamedToRunTimeSelectionTable(FieldValuePluginFunction, mqCellShapePluginFunction , name, mqCellShape);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mqCellShapePluginFunction::mqCellShapePluginFunction(
    const FieldValueExpressionDriver &parentDriver,
    const word &name
):
    FieldValuePluginFunction(
        parentDriver,
        name,
        word("volScalarField"),
        string("")
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void mqCellShapePluginFunction::doEvaluation()
{
    autoPtr<volScalarField> pShape(
        new volScalarField(
            IOobject(
                "cellShape",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar(dimless, 0),
            "zeroGradient"
        )
    );

    volScalarField &shape=pShape();

    // Construct shape recognizers
    hexMatcher hex;
    prismMatcher prism;
    wedgeMatcher wedge;
    pyrMatcher pyr;
    tetWedgeMatcher tetWedge;
    tetMatcher tet;

    forAll(shape,cellI) {
        if (hex.isA(mesh(),cellI)) {
            shape[cellI]=1;
        } else if (prism.isA(mesh(),cellI)) {
            shape[cellI]=2;
        } else if (wedge.isA(mesh(),cellI)) {
            shape[cellI]=3;
        } else if (pyr.isA(mesh(),cellI)) {
            shape[cellI]=4;
        } else if (tetWedge.isA(mesh(),cellI)) {
            shape[cellI]=5;
        } else if (tet.isA(mesh(),cellI)) {
            shape[cellI]=6;
        } else {
            shape[cellI]=0;
        }
    }

    shape.correctBoundaryConditions();

    result().setObjectResult(pShape);
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace

// ************************************************************************* //
