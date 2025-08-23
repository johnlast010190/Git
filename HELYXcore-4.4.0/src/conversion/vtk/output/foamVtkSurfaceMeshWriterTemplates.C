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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2016-2018 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::vtk::surfaceMeshWriter::getFaceField
(
    const SurfaceField<Type>& sfld
) const
{
    const polyBoundaryMesh& patches = sfld.mesh().boundaryMesh();

    const labelList& faceAddr = this->patch().addressing();

    auto tfld = tmp<Field<Type>>::New(faceAddr.size());
    auto iter = tfld.ref().begin();

    for (const label facei : faceAddr)
    {
        const label patchi = patches.whichPatch(facei);

        if (patchi == -1)
        {
            *iter = sfld[facei];
        }
        else
        {
            const label localFacei = facei - patches[patchi].start();
            *iter = sfld.boundaryField()[patchi][localFacei];
        }

        ++iter;
    }

    return tfld;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::vtk::surfaceMeshWriter::write
(
    const SurfaceField<Type>& field
)
{
    if (notState(outputState::CELL_DATA))
    {
        FatalErrorInFunction
            << "Bad writer state (" << stateNames[state_]
            << ") - should be (" << stateNames[outputState::CELL_DATA]
            << ") for field " << field.name() << nl << endl
            << exit(FatalError);
    }

    this->indirectPatchWriter::write(field.name(), getFaceField(field)());
}


template<class Type>
void Foam::vtk::surfaceMeshWriter::write
(
    const GeometricField<Type, faPatchField, areaMesh>& field
)
{
    if (notState(outputState::CELL_DATA))
    {
        FatalErrorInFunction
            << "Bad writer state (" << stateNames[state_]
            << ") - should be (" << stateNames[outputState::CELL_DATA]
            << ") for field " << field.name() << nl << endl
            << exit(FatalError);
    }

    this->indirectPatchWriter::write(field.name(), field.primitiveField());
}


// ************************************************************************* //
