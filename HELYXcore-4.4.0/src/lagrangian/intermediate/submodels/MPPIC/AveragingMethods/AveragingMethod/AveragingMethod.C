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
    (c) 2013-2017 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "submodels/MPPIC/AveragingMethods/AveragingMethod/AveragingMethod.H"
#include "meshes/pointMesh/pointMesh.H"
#include "fields/GeometricFields/pointFields/pointFieldsFwd.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::AveragingMethod<Type>::updateGrad()
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::AveragingMethod<Type>::AveragingMethod
(
    const IOobject& io,
    const dictionary& dict,
    const fvMesh& mesh,
    const labelList& size
)
:
    regIOobject(io),
    FieldField<Field, Type>(),
    dict_(dict),
    mesh_(mesh)
{
    forAll(size, i)
    {
        FieldField<Field, Type>::append
        (
            new Field<Type>(size[i], Zero)
        );
    }
}


template<class Type>
Foam::AveragingMethod<Type>::AveragingMethod
(
    const AveragingMethod<Type>& am
)
:
    regIOobject(am),
    FieldField<Field, Type>(am),
    dict_(am.dict_),
    mesh_(am.mesh_)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::AveragingMethod<Type>>
Foam::AveragingMethod<Type>::New
(
    const IOobject& io,
    const dictionary& dict,
    const fvMesh& mesh
)
{
    word averageType(dict.lookup(typeName));

    //Info<< "Selecting averaging method "
    //    << averageType << endl;

    const auto ctor =
        ctorTableLookup
        (
            "averaging method",
            dictionaryConstructorTable_(),
            averageType
        );
    return autoPtr<AveragingMethod<Type>>(ctor(io, dict, mesh));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::AveragingMethod<Type>::~AveragingMethod()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::AveragingMethod<Type>::average()
{
    updateGrad();
}


template<class Type>
void Foam::AveragingMethod<Type>::average
(
    const AveragingMethod<scalar>& weight
)
{
    updateGrad();

    *this /= max(weight, SMALL);
}


template<class Type>
bool Foam::AveragingMethod<Type>::writeData(Ostream& os) const
{
    return os.good();
}


template<class Type>
bool Foam::AveragingMethod<Type>::write(const bool write) const
{
    const pointMesh pointMesh_(mesh_);

    // point volumes
    Field<scalar> pointVolume(mesh_.nPoints(), 0);

    // output fields
    VolField<Type> cellValue
    (
        IOobject
        (
            this->name() + ":cellValue",
            this->time().timeName(),
            mesh_
        ),
        mesh_,
        dimensioned<Type>("zero", dimless, Zero)
    );
    VolField<TypeGrad> cellGrad
    (
        IOobject
        (
            this->name() + ":cellGrad",
            this->time().timeName(),
            mesh_
        ),
        mesh_,
        dimensioned<TypeGrad>("zero", dimless, Zero)
    );
    PointField<Type> pointValue
    (
        IOobject
        (
            this->name() + ":pointValue",
            this->time().timeName(),
            mesh_
        ),
        pointMesh_,
        dimensioned<Type>("zero", dimless, Zero)
    );
    PointField<TypeGrad> pointGrad
    (
        IOobject
        (
            this->name() + ":pointGrad",
            this->time().timeName(),
            mesh_
        ),
        pointMesh_,
        dimensioned<TypeGrad>("zero", dimless, Zero)
    );

    // tet-volume weighted sums
    forAll(mesh_.C(), celli)
    {
        const List<tetIndices> cellTets =
            polyMeshTetDecomposition::cellTetIndices(mesh_, celli);

        forAll(cellTets, tetI)
        {
            const tetIndices& tetIs = cellTets[tetI];
            const scalar v = tetIs.tet(mesh_).mag();

            cellValue[celli] += v*interpolate(mesh_.C()[celli], tetIs);
            cellGrad[celli] += v*interpolateGrad(mesh_.C()[celli], tetIs);

            const face& f = mesh_.faces()[tetIs.face()];
            labelList vertices(3);
            vertices[0] = f[tetIs.faceBasePt()];
            vertices[1] = f[tetIs.facePtA()];
            vertices[2] = f[tetIs.facePtB()];

            forAll(vertices, vertexI)
            {
                const label pointi = vertices[vertexI];

                pointVolume[pointi] += v;
                pointValue[pointi] +=
                    v*interpolate(mesh_.points()[pointi], tetIs);
                pointGrad[pointi] +=
                    v*interpolateGrad(mesh_.points()[pointi], tetIs);
            }
        }
    }

    // average
    cellValue.primitiveFieldRef() /= mesh_.V();
    cellGrad.primitiveFieldRef() /= mesh_.V();
    pointValue.primitiveFieldRef() /= pointVolume;
    pointGrad.primitiveFieldRef() /= pointVolume;

    // write
    if (!cellValue.write(write)) return false;
    if (!cellGrad.write(write)) return false;
    if (!pointValue.write(write)) return false;
    if (!pointGrad.write(write)) return false;

    return true;
}


// ************************************************************************* //
