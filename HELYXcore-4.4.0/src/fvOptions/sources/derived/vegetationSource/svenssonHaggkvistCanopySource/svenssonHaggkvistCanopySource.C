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
    (c) 2015-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "svenssonHaggkvistCanopySource.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "cfdTools/general/include/fvCFD.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "cfdTools/general/fvOptions/fvOption.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(svenssonHaggkvistCanopySource, 0);
    addToRunTimeSelectionTable
    (
        option,
        svenssonHaggkvistCanopySource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::svenssonHaggkvistCanopySource::svenssonHaggkvistCanopySource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    vegetationSource(name, modelType, dict, obr),
    CpEps1_(dimensionedScalar::lookupOrAddToDict("CpEps1", coeffs_, 1.8))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::svenssonHaggkvistCanopySource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (eqn.psi().name() == word(UName_))
    {
        const volVectorField& U = eqn.psi();
        eqn -=  fvm::Sp(canopy_()*mag(U), U);
    }
}


void Foam::fv::svenssonHaggkvistCanopySource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (eqn.psi().name() == word(UName_))
    {
        const volVectorField& U = eqn.psi();
        eqn -=  fvm::Sp(rho*canopy_() * mag(U), U);
    }
 }


void Foam::fv::svenssonHaggkvistCanopySource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (eqn.psi().name() == word("k"))
    {
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
        eqn += canopy_()*pow(mag(U),3);
    }
    else if (eqn.psi().name() == word("epsilon"))
    {
        const volScalarField& epsilon = eqn.psi();
        const volScalarField& k = obr_.lookupObject<volScalarField>("k");
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
        eqn += fvm::Sp(canopy_()*CpEps1_/k*pow(mag(U),3), epsilon);
    }
    else if (eqn.psi().name() == word("omega"))
    {
        const volScalarField& omega = eqn.psi();
        const volScalarField& k = obr_.lookupObject<volScalarField>("k");
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
        eqn += fvm::Sp(canopy_()*CpEps1_/k*pow(mag(U),3), omega);
    }
    else if (eqn.psi().name() == word(TName_))
    {
        const scalarField& V = mesh_.V();
        const volScalarField rho(getRho());
        const volScalarField& Cp = CpRef();

        scalarField& TSource = eqn.source();

        forAll(cells_, i)
        {
            const label celli = cells_[i];
            TSource[celli] -=LAD[i]*Qs[i]/(rho[celli]*Cp[celli])*V[celli];
        }
    }
    else if (eqn.psi().name() == word("w"))
    {
        const scalarField& V = mesh_.V();
        const volScalarField rho(getRho());

        scalarField& humiditySource = eqn.source();

        forAll(cells_, i)
        {
            const label celli = cells_[i];
            humiditySource[celli] -=
                LAD[i]*Ql[i]/(rho[celli]*lambdaVap_)*V[celli];
        }
    }
}


void Foam::fv::svenssonHaggkvistCanopySource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (eqn.psi().name() == word("k"))
    {
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
        eqn += rho*canopy_()*pow(mag(U), 3);
    }
    else if (eqn.psi().name() == word("epsilon"))
    {
        const volScalarField& epsilon = eqn.psi();
        const volScalarField& k = obr_.lookupObject<volScalarField>("k");
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
        eqn += fvm::Sp(rho*canopy_()*CpEps1_/k*pow(mag(U),3), epsilon);
    }
    else if (eqn.psi().name() == word("omega"))
    {
        const volScalarField& omega = eqn.psi();
        const volScalarField& k = obr_.lookupObject<volScalarField>("k");
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
        eqn += fvm::Sp(rho*canopy_()*CpEps1_/k*pow(mag(U),3), omega);
    }
    else if (eqn.psi().name() == heName_ || eqn.psi().name() == TName_)
    {
        const scalarField& V = mesh_.V();
        scalarField& heSource = eqn.source();
        forAll(cells_, i)
        {
            const label celli = cells_[i];
            heSource[celli] -=LAD[i]*Qs[i]*V[celli];
        }
    }
    else if (eqn.psi().name() == word("w"))
    {
        const scalarField& V = mesh_.V();
        scalarField& humiditySource = eqn.source();

        forAll(cells_, i)
        {
            const label celli = cells_[i];
            humiditySource[celli] -=LAD[i]*Ql[i]/lambdaVap_*V[celli];
        }
    }
}


// ************************************************************************* //
