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
    (c) 2014-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "segregated.H"
#include "../../../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianDragModels
{
    defineTypeNameAndDebug(segregated, 0);
    addToRunTimeSelectionTable(eulerianDragModel, segregated, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianDragModels::segregated::segregated
(
    const dictionary& dict,
    const eulerianPhasePair& pair,
    const bool registerObject
)
:
    eulerianDragModel(dict, pair, registerObject),
    m_("m", dimless, dict),
    n_("n", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianDragModels::segregated::~segregated()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::eulerianDragModels::segregated::CdRe() const
{
    FatalErrorInFunction
        << "Not implemented."
        << "Drag coefficient not defined for the segregated model."
        << exit(FatalError);

    return pair_.phase1().volFrac();
}


Foam::tmp<Foam::volScalarField> Foam::eulerianDragModels::segregated::K() const
{
    const fvMesh& mesh(pair_.phase1().mesh());

    const volScalarField& alpha1(pair_.phase1());
    const volScalarField& alpha2(pair_.phase2());

    tmp<volScalarField> trho1(pair_.phase1().rho());
    tmp<volScalarField> trho2(pair_.phase2().rho());
    const volScalarField& rho1 = trho1();
    const volScalarField& rho2 = trho2();

    tmp<volScalarField> tnu1(pair_.phase1().nu());
    tmp<volScalarField> tnu2(pair_.phase2().nu());

    const volScalarField& nu1(tnu1());
    const volScalarField& nu2(tnu2());

    volScalarField L
    (
        IOobject
        (
            "L",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimLength, 0),
        zeroGradientFvPatchField<scalar>::typeName
    );
    L.primitiveFieldRef() = cbrt(mesh.V());
    L.correctBoundaryConditions();

    volScalarField I
    (
        alpha1
       /max
        (
            alpha1 + alpha2,
            pair_.phase1().residualAlpha() + pair_.phase2().residualAlpha()
        )
    );
    volScalarField magGradI
    (
        max
        (
            mag(fvc::grad(I)),
            (pair_.phase1().residualAlpha() + pair_.phase2().residualAlpha())/L
        )
    );

    volScalarField muI
    (
        rho1*nu1*rho2*nu2
       /(rho1*nu1 + rho2*nu2)
    );
    volScalarField muAlphaI
    (
        alpha1*rho1*nu1*alpha2*rho2*nu2
       /(
           max(alpha1, pair_.phase1().residualAlpha())*rho1*nu1
         + max(alpha2, pair_.phase2().residualAlpha())*rho2*nu2
        )
    );

    volScalarField ReI
    (
        pair_.rho()
       *pair_.magUr()
       /(magGradI*muI)
    );

    volScalarField lambda(m_*ReI + n_*muAlphaI/muI);

    return lambda*sqr(magGradI)*muI;
}


Foam::tmp<Foam::surfaceScalarField> Foam::eulerianDragModels::segregated::Kf() const
{
    return fvc::interpolate(K());
}


// ************************************************************************* //
