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

#include "waveTheory.H"
#include "multiphaseThermo/multiphaseThermo.H"
#include "equationOfState/rhoConstLaw/rhoConstLaw.H"
#include "basicThermo/basicThermo.H"
#include "helyxSolve.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(waveTheory, 0);
defineRunTimeSelectionTable(waveTheory, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


waveTheory::waveTheory
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    IOdictionary
    (
        mesh_.thisDb().lookupObject<IOobject>("waveProperties")
    ),

    seaLevel_(lookup<scalar>("seaLevel")),

    g_
    (
        uniformDimensionedVectorField
        (
            mesh_.thisDb().lookupObject<uniformDimensionedVectorField>("g")
        ).value()
    ),

    direction_( g_/mag(g_) ),

    coeffDict_(subDict(subDictName + "Coeffs")),

    PI_( M_PI ),

    wind_( lookupOrDefault<vector>( "wind", vector::zero ) )
{
    const helyxSolve* helyxSolvePtr =
        mesh_.time().lookupObjectPtr<helyxSolve>(helyxSolve::typeName);

    if (helyxSolvePtr)
    {
        const multiphaseThermo* thermoPtr =
            mesh_.thisDb().lookupObjectPtr<multiphaseThermo>(basicThermo::dictName);

        if (thermoPtr)
        {
            // USF
            if ( !thermoPtr->thermos(waves2Foam::waterPhase()).isochoric())
            {
                FatalErrorInFunction
                    << "Only constant-density material properties are supported "
                    << "by waves2Foam. Please modify the equation of state to '"
                    << rhoConstLaw::typeName << "'."
                    << exit(FatalError);
            }
            rhoWater_ =
                gAverage(thermoPtr->thermos(waves2Foam::waterPhase()).rho()());
        }
        else
        {
            const basicThermo* thermoSfPtr = 
                &refCast<basicThermo>(basicThermo::lookupOrCreate(mesh_.thisDb()));

            if (!thermoSfPtr->isochoric())
            {
                FatalErrorInFunction
                    << "Only constant-density material properties are supported "
                    << "by waves2Foam. Please modify the equation of state to '"
                    << rhoConstLaw::typeName << "'."
                    << exit(FatalError);               
            }
            rhoWater_ = gAverage(thermoSfPtr->rho()());
        }
    }
    else
    {
        IOdictionary transProp
        (
            IOobject
            (
                "transportProperties",
                "constant",
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        dictionary sD(transProp.subDict(Foam::waves2Foam::waterPhase()));
        rhoWater_ = (dimensionedScalar(sD.lookup("rho"))).value();
    }
}


waveTheory::~waveTheory()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void waveTheory::checkWaveDirection(const vector& k) const
{
    if (Foam::mag(k & direction_) > SMALL)
    {
        FatalErrorIn("void waveTheory::checkWaveDirection(const vector& k)")
            << "The wave number " << k << " is not perpendicular to the \n"
            << "direction of the gravitational vector " << g_ << "\n"
            << endl << endl << exit(FatalError);
    }
}


void waveTheory::checkWaveDirection(const vectorField& k) const
{
    forAll(k, ki)
    {
        checkWaveDirection(k[ki]);
    }
}


scalarField waveTheory::eta
(
    const pointField& x,
    const scalar& time
) const
{
    scalarField temp(x.size(),0.0);

    forAll(x,ii)
    {
        temp[ii] = eta(x[ii],time);
    }

    return temp;
}


scalarField waveTheory::ddxPd
(
    const pointField& x,
    const scalar& time,
    const vectorField& unitVector
) const
{
    scalarField temp(x.size(),0.0);

    forAll(x,ii)
    {
        temp[ii] = ddxPd(x[ii],time, unitVector[ii]);
    }

    return temp;
}


vectorField waveTheory::U
(
    const pointField& x,
    const scalar& time
) const
{
    vectorField temp(x.size(),vector::zero);

    forAll(x,ii)
    {
        temp[ii] = U(x[ii],time);
    }

    return temp;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
