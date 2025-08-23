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
    (c) 2016 OpenCFD Ltd
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "pyrolysisModels/thermo/pyrolysisThermo.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "absorptionEmissionModels/absorptionEmissionModel/absorptionEmissionModel.H"
#include "finiteVolume/fvm/fvm.H"
#include "finiteVolume/fvc/fvcLaplacian.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
    defineTypeNameAndDebug(pyrolysisThermo, 0);
    addToRunTimeSelectionTable(pyrolysisModel, pyrolysisThermo, mesh);
    addToRunTimeSelectionTable(pyrolysisModel, pyrolysisThermo, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionModels::pyrolysisThermo::pyrolysisThermo
(
    const word& modelType,
    const fvMesh& mesh,
    const word& regionType
)
:
    pyrolysisModel(modelType, mesh, regionType)
{}


Foam::regionModels::pyrolysisThermo::pyrolysisThermo
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& regionType
)
:
    pyrolysisModel(modelType, mesh, dict, regionType)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionModels::pyrolysisThermo::~pyrolysisThermo()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::regionModels::pyrolysisThermo::evolveRegion()
{
    if (active())
    {
        Info<< "\nEvolving energy in region: " << regionMesh().name() << endl;
        volScalarField& he = (*solidThermo_).he();

        for (int nonOrth=0; nonOrth<=pimple_.correctNonOrthogonal(); nonOrth++)
        {
            tmp<volScalarField> talpha
            (
                (*solidThermo_).kappa()/(*solidThermo_).Cp()
            );
            volScalarField& alpha = talpha.ref();
            alpha.rename(IOobject::groupName("alpha", he.group()));

            fvScalarMatrix heEqn
            (
                fvm::ddt(rho(), he)
              - fvm::laplacian(alpha, he) + fvc::laplacian(alpha, he)
              - fvc::laplacian((*solidThermo_).kappa(), (*solidThermo_).T())
            );
            heEqn.relax();
            heEqn.solve();
        }
        (*solidThermo_).correct();
    }
}


// ************************************************************************* //
