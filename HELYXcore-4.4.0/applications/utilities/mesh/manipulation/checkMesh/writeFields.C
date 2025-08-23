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
    (c) 2014-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "writeFields.H"
#include "fields/fvPatchFields/basic/calculated/calculatedFvPatchFields.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void maxFaceToCell(const scalarField& faceData, volScalarField& cellData)
{
    const cellList& cells = cellData.mesh().cells();

    scalarField& cellFld = cellData.ref();

    cellFld = -GREAT;
    forAll(cells, cellI)
    {
        const cell& cFaces = cells[cellI];
        forAll(cFaces, i)
        {
            cellFld[cellI] = max(cellFld[cellI], faceData[cFaces[i]]);
        }
    }

    forAll(cellData.boundaryField(), patchI)
    {
        fvPatchScalarField& fvp = cellData.boundaryFieldRef()[patchI];

        fvp = fvp.patch().patchSlice(faceData);
    }
    cellData.correctBoundaryConditions();
}


void minFaceToCell(const scalarField& faceData, volScalarField& cellData)
{
    const cellList& cells = cellData.mesh().cells();

    scalarField& cellFld = cellData.ref();

    cellFld = GREAT;
    forAll(cells, cellI)
    {
        const cell& cFaces = cells[cellI];
        forAll(cFaces, i)
        {
            cellFld[cellI] = min(cellFld[cellI], faceData[cFaces[i]]);
        }
    }

    forAll(cellData.boundaryField(), patchI)
    {
        fvPatchScalarField& fvp = cellData.boundaryFieldRef()[patchI];

        fvp = fvp.patch().patchSlice(faceData);
    }
    cellData.correctBoundaryConditions();
}


void Foam::writeFields(const fvMesh& mesh, const HashSet<word>& selectedFields)
{
    if (selectedFields.empty())
    {
        return;
    }

    Info<< "Writing fields " << endl;

    // Cell type
    if (selectedFields.found("cellShapes"))
    {
        volScalarField shape
        (
            IOobject
            (
                "cellShapes",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar(dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        );
        const cellShapeList& cellShapes = mesh.cellShapes();
        forAll(cellShapes, cellI)
        {
            const cellModel& model = cellShapes[cellI].model();
            shape[cellI] = model.index();
        }
        shape.correctBoundaryConditions();
        Info<< "    Writing cell shape (hex, tet etc.) to " << shape.name()
            << endl;
        shape.write();
    }

    if (selectedFields.found("cellVolume"))
    {
        volScalarField V
        (
            IOobject
            (
                "cellVolume",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar(dimVolume, 0),
            calculatedFvPatchScalarField::typeName
        );
        V.ref() = mesh.V();
        Info<< "    Writing cell volume to " << V.name() << endl;
        V.write();
    }

    if (selectedFields.found("cellCentres"))
    {
        volVectorField C
        (
            IOobject
            (
                "cellCentres",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh,
            dimensionedVector("cellCentres", dimLength, vector::zero),
            calculatedFvPatchVectorField::typeName
        );
        C.ref() = mesh.C();
        Info<< "    Writing cell centres to " << C.name() << endl;
        C.write();
    }

    Info<< endl;
}


// ************************************************************************* //
