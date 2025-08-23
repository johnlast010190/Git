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
    (c) held by original author

\*---------------------------------------------------------------------------*/

#include "relaxationZone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(relaxationZone, 0);

// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * //


relaxationZone::relaxationZone
(
    const fvMesh& mesh,
    volVectorField& U,
    volScalarField& alpha
)
:
    mesh_(mesh),
    U_(U),
    alpha_(alpha),

    relaxNames_((mesh_.thisDb().lookupObject<IOdictionary>("waveProperties"))
        .lookup("relaxationNames")),

    relaxSchemePtr_(relaxNames_.size()),

    relaxationWeightsMomentum_(nullptr),

    relaxationWeightsSource_(nullptr),

    targetAlpha_(nullptr),

    targetVelocity_(nullptr)


{
    forAll(relaxNames_, relaxi)
    {
        relaxSchemePtr_[relaxi] = relaxationSchemes::relaxationScheme::
            New(relaxNames_[relaxi], mesh_, U_, alpha_);
    }
}


// * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * //


void relaxationZone::resetTargetFields()
{
    if (relaxationWeightsMomentum_ != nullptr)
    {
#if EXTBRANCH==1
        (*relaxationWeightsMomentum_).internalField() = 0.0;
#elif OFPLUSBRANCH==1
    #if OFVERSION<1706
        (*relaxationWeightsMomentum_).internalField() = 0.0;
    #else
        (*relaxationWeightsMomentum_).primitiveFieldRef() = 0.0;
    #endif
#else
    #if OFVERSION<400
        (*relaxationWeightsMomentum_).internalField() = 0.0;
    #else
        (*relaxationWeightsMomentum_).ref() = 0.0;
    #endif
#endif
    }

    if (relaxationWeightsSource_ != nullptr)
    {
#if EXTBRANCH==1
        (*relaxationWeightsSource_).internalField() = 1.0;
#elif OFPLUSBRANCH==1
    #if OFVERSION<1706
        (*relaxationWeightsSource_).internalField() = 1.0;
    #else
        (*relaxationWeightsSource_).primitiveFieldRef() = 1.0;
    #endif
#else
    #if OFVERSION<400
        (*relaxationWeightsSource_).internalField() = 1.0;
    #else
        (*relaxationWeightsSource_).ref() = 1.0;
    #endif
#endif
    }

    if (targetAlpha_ != nullptr)
    {
#if EXTBRANCH==1
        (*targetAlpha_).internalField() = 0.0;
#elif OFPLUSBRANCH==1
    #if OFVERSION<1706
        (*targetAlpha_).internalField() = 0.0;
    #else
        (*targetAlpha_).primitiveFieldRef() = 0.0;
    #endif
#else
    #if OFVERSION<400
        (*targetAlpha_).internalField() = 0.0;
    #else
        (*targetAlpha_).ref() = 0.0;
    #endif
#endif
    }

    if (targetVelocity_ != nullptr)
    {
#if EXTBRANCH==1
        (*targetVelocity_).internalField() = vector::zero;
#elif OFPLUSBRANCH==1
    #if OFVERSION<1706
        (*targetVelocity_).internalField() = vector::zero;
    #else
        (*targetVelocity_).primitiveFieldRef() = vector::zero;
    #endif
#else
    #if OFVERSION<400
        (*targetVelocity_).internalField() = vector::zero;
    #else
        (*targetVelocity_).ref() = vector::zero;
    #endif
#endif
    }
}


void relaxationZone::correctBoundaries()
{
    alpha_.correctBoundaryConditions();

    U_.correctBoundaryConditions();

    if (relaxationWeightsMomentum_ != nullptr)
    {
        (*relaxationWeightsMomentum_).correctBoundaryConditions();
    }

    if (relaxationWeightsSource_ != nullptr)
    {
        (*relaxationWeightsSource_).correctBoundaryConditions();
    }

    if (targetAlpha_ != nullptr)
    {
        (*targetAlpha_).correctBoundaryConditions();
    }

    if (targetVelocity_ != nullptr)
    {
        (*targetVelocity_).correctBoundaryConditions();
    }
}


void relaxationZone::correct()
{
    scalar preTime = mesh_.time().elapsedCpuTime();

    resetTargetFields();

    forAll(relaxSchemePtr_, relaxi)
    {
        relaxSchemePtr_[relaxi]->correct();
    }

    correctBoundaries();

    Info<< "Relaxing time: " << mesh_.time().elapsedCpuTime() - preTime
         << " s" << endl;
}


tmp<volScalarField> relaxationZone::numericalBeach()
{
    // Return field
    tmp<volScalarField> tartificialViscotity
    (
        volScalarField::New
        (
            "muArti",
            mesh_,
            dimensionedScalar(dimless/dimTime, 0),
            "fixedValue"
        )
    );

#if EXTBRANCH==1
    volScalarField& artificialViscosity(tartificialViscotity());
#elif OFPLUSBRANCH==1
    #if OFVERSION<1606
        volScalarField& artificialViscosity(tartificialViscotity());
    #else
        volScalarField& artificialViscosity(tartificialViscotity.ref());
    #endif
#else
    #if OFVERSION<400
        volScalarField& artificialViscosity(tartificialViscotity());
    #else
        volScalarField& artificialViscosity(tartificialViscotity.ref());
    #endif
#endif

    forAll(relaxSchemePtr_, relaxi)
    {
        relaxSchemePtr_[relaxi]->numericalBeach( artificialViscosity );
    }

    return tartificialViscotity;
}


const volScalarField& relaxationZone::relaxationWeightsMomentum() const
{
    if (!relaxationWeightsMomentum_)
    {
        relaxationWeightsMomentum_ = new volScalarField
        (
            IOobject
            (
                "relaxationWeightsMomentum",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, 0),
            "zeroGradient"
        );
    }

    return *relaxationWeightsMomentum_;
}


const volScalarField& relaxationZone::relaxationWeightsSource() const
{
    if (!relaxationWeightsSource_)
    {
        relaxationWeightsSource_ = new volScalarField
        (
            IOobject
            (
                "relaxationWeightsSource",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, 1),
            "zeroGradient"
        );
    }

    return *relaxationWeightsSource_;
}


const volScalarField& relaxationZone::targetAlphaField() const
{
    if (!targetAlpha_)
    {
        targetAlpha_ = new volScalarField
        (
            IOobject
            (
                "targetAlphaField",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, 0),
            "zeroGradient"
        );
    }

    return *targetAlpha_;
}


const volVectorField& relaxationZone::targetVelocityField() const
{
    if (!targetVelocity_)
    {
        targetVelocity_ = new volVectorField
        (
            IOobject
            (
                "targetVelocityField",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedVector(dimVelocity, vector::zero),
            "zeroGradient"
        );
    }

    return *targetVelocity_;
}


} // end namespace Foam
