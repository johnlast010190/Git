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
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "surfaceHasSurfacePluginFunction.H"
#include "FieldValueExpressionDriver.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

namespace Foam {

defineTypeNameAndDebug(surfaceHasSurfacePluginFunction,0);
addNamedToRunTimeSelectionTable(FieldValuePluginFunction, surfaceHasSurfacePluginFunction , name, surfaceHasSurface);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

surfaceHasSurfacePluginFunction::surfaceHasSurfacePluginFunction(
    const FieldValueExpressionDriver &parentDriver,
    const word &name
):
    GeneralSurfacesPluginFunction(
        parentDriver,
        name,
        "volLogicalField",
        string("surfaceName primitive word")
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void surfaceHasSurfacePluginFunction::doEvaluation()
{
    autoPtr<volScalarField> pHasSurface(
        new volScalarField(
            IOobject(
                "surfaceHasSurfaceInCell",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar(dimless, 0)
        )
    );

    const labelList &cells=meshCells();

    forAll(cells,i) {
        const label cellI=cells[i];

        if (cellI>=0) {
            pHasSurface()[cellI]=1;
        }
    }

    pHasSurface->correctBoundaryConditions();

    result().setObjectResult(pHasSurface);
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace

// ************************************************************************* //
