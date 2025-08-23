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
    (c) 2017 OpenCFD Ltd.
    (c) 2022-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/fv/fv.H"
#include "containers/HashTables/HashTable/HashTable.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "fvMesh/fvPatches/constraint/cyclicAMI/cyclicAMIFvPatch.H"
#include "fields/fvPatchFields/basic/blended/blendedFvPatchFields.H"
#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
tmp<ddtScheme<Type>> ddtScheme<Type>::New
(
    const fvMesh& mesh,
    Istream& schemeData
)
{
    if (fv::debug)
    {
        InfoInFunction << "Constructing ddtScheme<Type>" << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Ddt scheme not specified" << endl << endl
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    const auto ctor =
        ctorTableLookup
        (
            "ddt scheme",
            IstreamConstructorTable_(),
            schemeName
        );
    return ctor(mesh, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
ddtScheme<Type>::~ddtScheme()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolField<Type>> ddtScheme<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    NotImplemented;
}


template<class Type>
tmp<fvMatrix<Type>> ddtScheme<Type>::fvmDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    NotImplemented;
}


template<class Type>
tmp<SurfaceField<Type>> ddtScheme<Type>::fvcDdt
(
    const SurfaceField<Type>& sf
)
{
    NotImplemented;
}


template<class Type>
scalar ddtScheme<Type>::getDdtPhiCoeff(const fvMesh& mesh)
{
    if (!ddtPhiCoeff_.valid())
    {
        const solutionControl* solControl =
            mesh.lookupObjectPtr<solutionControl>(solutionControl::typeName);
        ddtPhiCoeff_.reset
        (
            solControl && solControl->dict().found("ddtPhiCoeff")
          ? Function1<scalar>::New("ddtPhiCoeff", solControl->dict()).ptr()
          : new Function1Types::Constant<scalar>("ddtPhiCoeff", scalar(-1))
        );
    }

    const scalar ddtCoeff = ddtPhiCoeff_->value(mesh.time().value());
    DebugInFunction << "Setting ddtPhiCoeff to " << ddtCoeff << endl;

    return ddtCoeff;
}


template<class Type>
tmp<surfaceScalarField> ddtScheme<Type>::fvcDdtPhiCoeff
(
    const VolField<Type>& U,
    const fluxFieldType& phi,
    const fluxFieldType& phiCorr
)
{
    tmp<surfaceScalarField> tddtCouplingCoeff
    (
        surfaceScalarField::New
        (
            "ddtCouplingCoeff",
            U.mesh(),
            dimensionedScalar(dimless, 1)
        )
    );

    surfaceScalarField& ddtCouplingCoeff = tddtCouplingCoeff.ref();

    scalar ddtPhiCoeff(getDdtPhiCoeff(U.mesh()));

    if (ddtPhiCoeff < 0)
    {
        ddtCouplingCoeff -= min
        (
            mag(phiCorr)
           /(mag(phi) + dimensionedScalar(phi.dimensions(), SMALL)),
            scalar(-ddtPhiCoeff)
        );
    }
    else
    {
        ddtCouplingCoeff =
            dimensionedScalar("ddtPhiCoeff", dimless, ddtPhiCoeff);
    }

    surfaceScalarField::Boundary& ccbf = ddtCouplingCoeff.boundaryFieldRef();

    forAll(U.boundaryField(), patchi)
    {
        if (isA<blendedFvPatchVectorField>(U.boundaryField()[patchi]))
        {
            List<bool> fixesValues(U.boundaryField()[patchi].fixesValues());

            forAll(ccbf[patchi], facei)
            {
                if (fixesValues[facei])
                {
                    ccbf[patchi][facei] = 0.0;
                }
            }
        }
        else if
        (
            U.boundaryField()[patchi].fixesValue()
         || isA<cyclicAMIFvPatch>(mesh().boundary()[patchi])
        )
        {
            ccbf[patchi] = 0.0;
        }
    }

    if (debug > 1)
    {
        InfoInFunction
            << "ddtCouplingCoeff mean max min = "
            << gAverage(ddtCouplingCoeff.primitiveField())
            << " " << gMax(ddtCouplingCoeff.primitiveField())
            << " " << gMin(ddtCouplingCoeff.primitiveField())
            << endl;
    }

    return tddtCouplingCoeff;
}


template<class Type>
tmp<surfaceScalarField> ddtScheme<Type>::fvcDdtPhiCoeff
(
    const VolField<Type>& U,
    const fluxFieldType& phi
)
{
    return fvcDdtPhiCoeff(U, phi, phi - fvc::dotInterpolate(mesh().Sf(), U));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
