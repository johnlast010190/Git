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
    (c) 2016-2023 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "VOFState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "multiphaseThermo/multiphaseThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace stateFunctions
{

defineTypeNameAndDebug(VOFState, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

VOFState::VOFState
(
    word region,
    const dictionary& input,
    const dictionary& defaults,
    const stateFunction& master,
    const stateIndex& index,
    word meshName
)
:
    multiphaseMulticomponentState
    (
        region,
        input,
        defaults,
        master,
        index,
        meshName
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void VOFState::initialise()
{
    // User's explicitly defined passive phase, or first phase in two-phase case
    word passivePhaseName =
        input().subOrEmptyDict("materialProperties").lookupOrDefault
        (
            "passivePhase",
            phaseNames().size() == 2 ? phaseNames()[1] : word::null
        );

    forAll(phaseNames(), pi)
    {
        // Don't create a fieldMap for the passive phase
        if (phaseNames()[pi] == passivePhaseName) continue;

        // fieldMaps
        stateDict_.subDict("fieldMaps").add
        (
            word("alpha." + phaseNames()[pi]),
            input().subDict("fieldMaps").lookup("alpha")
        );
    }
    // Remove old one from input to prevent it being merged later
    input().subDict("fieldMaps").remove("alpha");

    multiphaseMulticomponentState::initialise();

    // Check that we have a multiphase material model
    dictionary& matDict = system().subDict("materialProperties");
    const word matType(matDict.lookup("materialType"));
    if (matType != "multiphase")
    {
        FatalErrorInFunction
            << "The 'VOF' state requires a multiphase material"
            << nl << exit(FatalError);
    }

    if (passivePhaseName != word::null)
    {
        matDict.add("passivePhase", passivePhaseName);
    }

    forAll(phaseNames(), pi)
    {
        // find grad(alpha) in defaults and write to state using proper phase name
        addAlphaTypeGrad(phaseNames()[pi]);
    }

    // Merge binary models from defaults
    mergeBinaryPhaseModels("surfaceTensionModel", false);
    mergeBinaryPhaseModels("phaseChangeModel", true, false);
}


void VOFState::finalise()
{
    multiphaseMulticomponentState::finalise();

    // clean-up
    if (system().found("fvSchemes"))
    {
        if (system().subDict("fvSchemes").found("gradSchemes"))
        {
            system().subDict("fvSchemes")
                .subDict("gradSchemes").remove("grad(alpha)");
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// ************************************************************************* //
