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
    (c) 2014-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "eulerianPhasePair.H"
#include "binaryPhaseModels/surfaceTensionModels/surfaceTensionModel/surfaceTensionModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::eulerianPhasePair::EoH
(
    const volScalarField& d
) const
{
    if (gPtr_)
    {
        return
            mag(dispersed().rho() - continuous().rho())*mag(*gPtr_)*sqr(d)
           /sigma();
    }
    else
    {
        return volScalarField::New
        (
            "EoH",
            d.db(),
            phase1().mesh(),
            dimensionedScalar(dimless, 0)
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianPhasePair::eulerianPhasePair
(
    const eulerianPhaseModel& phase1,
    const eulerianPhaseModel& phase2,
    const bool ordered
)
:
    phasePairKey(phase1.name(), phase2.name(), ordered),
    phase1_(phase1),
    phase2_(phase2),
    gPtr_(phase1.mesh().lookupObjectPtr<uniformDimensionedVectorField>("g"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianPhasePair::~eulerianPhasePair()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::eulerianPhaseModel& Foam::eulerianPhasePair::dispersed() const
{
    FatalErrorInFunction
        << "Requested dispersed phase from an unordered pair."
        << exit(FatalError);

    return phase1_;
}


const Foam::eulerianPhaseModel& Foam::eulerianPhasePair::continuous() const
{
    FatalErrorInFunction
        << "Requested continuous phase from an unordered pair."
        << exit(FatalError);

    return phase1_;
}


Foam::tmp<Foam::volScalarField> Foam::eulerianPhasePair::rho() const
{
    return
        phase1().volFrac()*phase1().rho() + phase2().volFrac()*phase2().rho();
}


Foam::tmp<Foam::volScalarField> Foam::eulerianPhasePair::magUr() const
{
    return mag(phase1().U() - phase2().U());
}


Foam::tmp<Foam::volVectorField> Foam::eulerianPhasePair::Ur() const
{
    return dispersed().U() - continuous().U();
}


Foam::tmp<Foam::volScalarField> Foam::eulerianPhasePair::Re() const
{
    return magUr()*dispersed().d()/continuous().nu();
}


Foam::tmp<Foam::volScalarField> Foam::eulerianPhasePair::Pr() const
{
    return
         continuous().nu()
        *continuous().thermo().Cpv()
        *continuous().rho()
        /continuous().kappa();
}


Foam::tmp<Foam::volScalarField> Foam::eulerianPhasePair::Eo() const
{
    return EoH(dispersed().d());
}


Foam::tmp<Foam::volScalarField> Foam::eulerianPhasePair::EoH1() const
{
    return EoH(dispersed().d()*cbrt(1 + 0.163*pow(Eo(), 0.757)));
}


Foam::tmp<Foam::volScalarField> Foam::eulerianPhasePair::EoH2() const
{
    return EoH(dispersed().d()/cbrt(E()));
}


Foam::tmp<Foam::volScalarField> Foam::eulerianPhasePair::sigma() const
{
    word sigmaName =
        IOobject::groupName
        (
            surfaceTensionModel::typeName,
            phasePairKey(first(), second(), false).name()
        );
    if (phase1().mesh().foundObject<surfaceTensionModel>(sigmaName))
    {
        return
            phase1().mesh().lookupObject<surfaceTensionModel>
            (
                sigmaName
            ).sigma();
    }
    else
    {
        return
            phase1().mesh().lookupObject<surfaceTensionModel>
            (
                IOobject::groupName
                (
                    surfaceTensionModel::typeName,
                    phasePairKey(second(), first(), false).name()
                )
            ).sigma();
    }
}


Foam::tmp<Foam::volScalarField> Foam::eulerianPhasePair::Mo() const
{
    if (gPtr_)
    {
        return
            mag(*gPtr_)
           *continuous().nu()
           *pow3
            (
                continuous().nu()
               *continuous().rho()
               /sigma()
            );
    }
    else
    {
        return volScalarField::New
        (
            "Mo",
            continuous().rho()().db(),
            phase1().mesh(),
            dimensionedScalar(dimless, 0)
        );
    }
}


Foam::tmp<Foam::volScalarField> Foam::eulerianPhasePair::Ta() const
{
    return Re()*pow(Mo(), 0.23);
}


Foam::tmp<Foam::volScalarField> Foam::eulerianPhasePair::E() const
{
    FatalErrorInFunction
        << "Requested aspect ratio of the dispersed phase in an unordered pair"
        << exit(FatalError);

    return phase1().volFrac();
}


// ************************************************************************* //
