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
    (c) 2016-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/constrainHbyA/constrainHbyA.H"
#include "fields/volFields/volFields.H"
#include "fields/fvPatchFields/derived/fixedFluxExtrapolatedPressure/fixedFluxExtrapolatedPressureFvPatchScalarField.H"
#include "fields/fvPatchFields/basic/blended/blendedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::constrainHbyA
(
    const tmp<volVectorField>& tHbyA,
    const volVectorField& U,
    const volScalarField& p
)
{
    tmp<volVectorField> tHbyANew;

    if (tHbyA.isTmp())
    {
        tHbyANew = tHbyA;
        tHbyANew.ref().rename(IOobject::groupName("HbyA", U.group()));
    }
    else
    {
        tHbyANew =
            volVectorField::New(IOobject::groupName("HbyA", U.group()), tHbyA);
    }

    volVectorField& HbyA = tHbyANew.ref();
    volVectorField::Boundary& HbyAbf = HbyA.boundaryFieldRef();

    forAll(U.boundaryField(), patchi)
    {
        if (isA<blendedFvPatchVectorField>(U.boundaryField()[patchi]))
        {
            const boolList assignables(U.boundaryField()[patchi].assignables());
            forAll(U.boundaryField()[patchi], facei)
            {
                if (!assignables[facei])
                {
                    HbyAbf[patchi][facei] = U.boundaryField()[patchi][facei];
                }
            }
        }
        else if
        (
           !U.boundaryField()[patchi].assignable()
        && !isA<fixedFluxExtrapolatedPressureFvPatchScalarField>
            (
                p.boundaryField()[patchi]
            )
        )
        {
            HbyAbf[patchi] = U.boundaryField()[patchi];
        }
    }

    return tHbyANew;
}


Foam::tmp<Foam::surfaceScalarField> Foam::constrainPhiHbyA
(
    const tmp<surfaceScalarField>& tphiHbyA,
    const volVectorField& U,
    const volScalarField& p
)
{
    tmp<surfaceScalarField> tphiHbyANew;

    if (tphiHbyA.isTmp())
    {
        tphiHbyANew = tphiHbyA;
        tphiHbyANew.ref().rename(IOobject::groupName("phiHbyA", U.group()));
    }
    else
    {
        tphiHbyANew = surfaceScalarField::New
        (
            IOobject::groupName("phiHbyA", U.group()),
            tphiHbyA
        );
    }

    surfaceScalarField& phiHbyA = tphiHbyANew.ref();
    surfaceScalarField::Boundary& phiHbyAbf = phiHbyA.boundaryFieldRef();

    forAll(U.boundaryField(), patchi)
    {
        if (isA<blendedFvPatchVectorField>(U.boundaryField()[patchi]))
        {
            const boolList assignables(U.boundaryField()[patchi].assignables());
            forAll(U.boundaryField()[patchi], facei)
            {
                if (!assignables[facei])
                {
                    phiHbyAbf[patchi][facei] =
                        U.mesh().Sf().boundaryField()[patchi][facei]
                      & U.boundaryField()[patchi][facei];
                }
            }
        }
        else if
        (
           !U.boundaryField()[patchi].assignable()
        && !isA<fixedFluxExtrapolatedPressureFvPatchScalarField>
            (
                p.boundaryField()[patchi]
            )
        )
        {
            phiHbyAbf[patchi] =
                U.mesh().Sf().boundaryField()[patchi]
              & U.boundaryField()[patchi];
        }
    }

    return tphiHbyANew;
}


// ************************************************************************* //
