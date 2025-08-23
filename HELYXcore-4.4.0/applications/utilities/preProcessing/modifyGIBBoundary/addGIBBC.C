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
    (c) 1991-2008 OpenCFD Ltd.
    (c) 2010-2016 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/boundarySetup/boundarySetup.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/directFvPatchFieldMapper.H"

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

template<class T>
void Foam::addGIBBC
(
    const dictionary& bf,
    const word fieldName,
    const fvMesh& initMesh,
    const fvMesh& updatedMesh
)
{
    typedef VolField<T> GeometricFieldType;

    bool found = initMesh.foundObject<GeometricFieldType>(fieldName);
    if (found)
    {
        GeometricFieldType& field =
            const_cast<GeometricFieldType&>
            (
                initMesh.lookupObject<GeometricFieldType>(fieldName)
            );

        PtrList<fvPatchField<T>> nubf(updatedMesh.boundary().size());

        //-- create existing new fvPatchFields
        //   Use this ::New function because the initMesh has changed and patch_
        //   reference inside the BCs is not correct
        forAll(field.boundaryField(), pI)
        {
            labelList map(identity(updatedMesh.boundary()[pI].size()));
            nubf.set
            (
                pI,
                fvPatchField<T>::New
                (
                    field.boundaryField()[pI],
                    initMesh.boundary()[pI],
                    field.internalField(),
                    directFvPatchFieldMapper(map)
                )
            );

        }

        //-- create GIB BCs
        int nI = field.boundaryField().size();
        forAllConstIter(IDLList<entry>, bf, iter)
        {
            nubf.set
            (
                nI
                ,
                    fvPatchField<T>::New
                    (
                        updatedMesh.boundary()[nI],
                        field.internalField(),
                        iter().dict()
                    )
            );
            nI++;
        }


        //-- create GIB BCs
        typename VolField<T>::
            Boundary newBCs
            (
                updatedMesh.boundary(),
                field.internalField(),
                nubf
            );

        GeometricFieldType fieldNewMesh(
            IOobject
            (
                field.name(),
                updatedMesh.time().timeName(),
                updatedMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            updatedMesh,
            field.dimensions(),
            field.internalField(),
            newBCs
        );
        fieldNewMesh.write();
    }
    else
    {
    }
}

// ************************************************************************* //
