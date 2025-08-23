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
    (c) 2016-2017 OpenCFD Ltd.
    (c) 2019-2023 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "jouleHeatingSource.H"
#include "solverObjects/solverOption/SolverOption.H"
#include "fvMatrices/fvMatrices.H"
#include "finiteVolume/fvm/fvmLaplacian.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchField.H"
#include "basicThermo/basicThermo.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(jouleHeatingSource, 0);
}
}

makeFvSolverOption(jouleHeatingSource);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::jouleHeatingSource::jouleHeatingSource
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::jouleHeatingSource::~jouleHeatingSource()
{}


// * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * //

bool Foam::fv::jouleHeatingSource::initialise()
{
    return true;
}


void Foam::fv::jouleHeatingSource::getSourceGraph
(
    wordList& fieldNames,
    HashTable<wordList>& sourceDependencies
)
{
    // Set the field name to that of the energy field from which the temperature
    // is obtained

    const basicThermo& thermo =
        obr_.lookupObject<basicThermo>(basicThermo::dictName);

    fieldNames.setSize(1, thermo.heT().name());

    sourceDependencies.insert(fieldNames[0], {"electrical_V"});
}


void Foam::fv::jouleHeatingSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    DebugInformation<< name() << ": applying source to " << eqn.psi().name() << endl;

    // Add the Joule heating source term

    const volScalarField& V =
        obr_.lookupObject<volScalarField>("electrical_V");
    const volVectorField gradV(fvc::grad(V));

    if (mesh().foundObject<volSymmTensorField>("electrical_sigma"))
    {
        // Anisotropic conductivity
        const volSymmTensorField& sigma =
            obr_.lookupObject<volSymmTensorField>("electrical_sigma");
        eqn += (sigma & gradV) & gradV;
    }
    else
    {
        const volScalarField& sigma =
            obr_.lookupObject<volScalarField>("electrical_sigma");
        eqn += (sigma*gradV) & gradV;
    }
}


// ************************************************************************* //
