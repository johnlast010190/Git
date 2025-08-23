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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "psiuMulticomponentThermo/psiuMulticomponentThermo.H"
#include "fvMesh/fvMesh.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "derivedFvPatchFields/fixedUnburntEnthalpy/fixedUnburntEnthalpyFvPatchScalarField.H"
#include "derivedFvPatchFields/gradientUnburntEnthalpy/gradientUnburntEnthalpyFvPatchScalarField.H"
#include "derivedFvPatchFields/mixedUnburntEnthalpy/mixedUnburntEnthalpyFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(psiuMulticomponentThermo, 0);
    defineRunTimeSelectionTable(psiuMulticomponentThermo, objectRegistry);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::psiuMulticomponentThermo::heuBoundaryTypes()
{
    const volScalarField::Boundary& tbf = this->Tu().boundaryField();

    wordList hbt = tbf.types();

    forAll(tbf, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = fixedUnburntEnthalpyFvPatchScalarField::typeName;
        }
        else if
        (
            isA<zeroGradientFvPatchScalarField>(tbf[patchi])
         || isA<fixedGradientFvPatchScalarField>(tbf[patchi])
        )
        {
            hbt[patchi] = gradientUnburntEnthalpyFvPatchScalarField::typeName;
        }
        else if (isA<mixedFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = mixedUnburntEnthalpyFvPatchScalarField::typeName;
        }
    }

    return hbt;
}


void Foam::psiuMulticomponentThermo::heuBoundaryCorrection(volScalarField& heu)
{
    volScalarField::Boundary& hbf = heu.boundaryFieldRef();

    forAll(hbf, patchi)
    {
        if (isA<gradientUnburntEnthalpyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<gradientUnburntEnthalpyFvPatchScalarField>(hbf[patchi])
                .gradient() = hbf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedUnburntEnthalpyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<mixedUnburntEnthalpyFvPatchScalarField>(hbf[patchi])
                .refGrad() = hbf[patchi].fvPatchField::snGrad();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::psiuMulticomponentThermo::implementation::implementation
(
    const objectRegistry& obr,
    const word& phaseName
)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::psiuMulticomponentThermo>
Foam::psiuMulticomponentThermo::New
(
    const objectRegistry& obr,
    const word& phaseName
)
{
    return basicThermo::New<psiuMulticomponentThermo>(obr, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::psiuMulticomponentThermo::~psiuMulticomponentThermo()
{}


Foam::psiuMulticomponentThermo::implementation::~implementation()
{}


// ************************************************************************* //
