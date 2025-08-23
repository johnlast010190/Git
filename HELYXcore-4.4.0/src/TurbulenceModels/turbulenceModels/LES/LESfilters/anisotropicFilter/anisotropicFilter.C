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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "anisotropicFilter.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "finiteVolume/fvc/fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(anisotropicFilter, 0);
    addToRunTimeSelectionTable(LESfilter, anisotropicFilter, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::anisotropicFilter::anisotropicFilter
(
    const fvMesh& mesh,
    scalar widthCoeff
)
:
    LESfilter(mesh),
    widthCoeff_(widthCoeff),
    coeff_
    (
        IOobject
        (
            "anisotropicFilterCoeff",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedVector(dimLength*dimLength, Zero),
        calculatedFvPatchVectorField::typeName
    )
{
    for (direction d=0; d<vector::nComponents; d++)
    {
        coeff_.primitiveFieldRef().replace
        (
            d,
            (1/widthCoeff_)*
            sqr
            (
                2.0*mesh.V()
               /fvc::surfaceSum(mag(mesh.Sf().component(d)))().primitiveField()
            )
        );
    }
}


Foam::anisotropicFilter::anisotropicFilter
(
    const fvMesh& mesh,
    const dictionary& bd
)
:
    LESfilter(mesh),
    widthCoeff_
    (
        bd.optionalSubDict(type() + "Coeffs").lookup<scalar>("widthCoeff")
    ),
    coeff_
    (
        IOobject
        (
            "anisotropicFilterCoeff",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedVector(dimLength*dimLength, Zero),
        calculatedFvPatchScalarField::typeName
    )
{
    for (direction d=0; d<vector::nComponents; d++)
    {
        coeff_.primitiveFieldRef().replace
        (
            d,
            (1/widthCoeff_)*
            sqr
            (
                2.0*mesh.V()
               /fvc::surfaceSum(mag(mesh.Sf().component(d)))().primitiveField()
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::anisotropicFilter::read(const dictionary& bd)
{
    widthCoeff_ =
        bd.optionalSubDict(type() + "Coeffs").lookup<scalar>("widthCoeff");
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::anisotropicFilter::operator()
(
    const tmp<volScalarField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volScalarField> tmpFilteredField =
        unFilteredField
      + (
           coeff_
         & fvc::surfaceIntegrate
           (
               mesh().Sf()
              *fvc::snGrad(unFilteredField())
           )
        );

    unFilteredField.clear();

    return tmpFilteredField;
}


Foam::tmp<Foam::volVectorField> Foam::anisotropicFilter::operator()
(
    const tmp<volVectorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volVectorField> tmpFilteredField =
        unFilteredField
      + (
           coeff_
         & fvc::surfaceIntegrate
           (
               mesh().Sf()
              *fvc::snGrad(unFilteredField())
           )
        );

    unFilteredField.clear();

    return tmpFilteredField;
}


Foam::tmp<Foam::volSymmTensorField> Foam::anisotropicFilter::operator()
(
    const tmp<volSymmTensorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volSymmTensorField> tmpFilteredField
    (
        volSymmTensorField::New
        (
            "anisotropicFilteredSymmTensorField",
            mesh(),
            unFilteredField().dimensions()
        )
    );

    for (direction d=0; d<symmTensor::nComponents; d++)
    {
        tmpFilteredField.ref().replace
        (
            d,
            anisotropicFilter::operator()(unFilteredField().component(d))
        );
    }

    unFilteredField.clear();

    return tmpFilteredField;
}


Foam::tmp<Foam::volTensorField> Foam::anisotropicFilter::operator()
(
    const tmp<volTensorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volTensorField> tmpFilteredField
    (
        volTensorField::New
        (
            "anisotropicFilteredTensorField",
            mesh(),
            unFilteredField().dimensions()
        )
    );

    for (direction d=0; d<tensor::nComponents; d++)
    {
        tmpFilteredField.ref().replace
        (
            d, anisotropicFilter::operator()(unFilteredField().component(d))
        );
    }

    unFilteredField.clear();

    return tmpFilteredField;
}


// ************************************************************************* //
