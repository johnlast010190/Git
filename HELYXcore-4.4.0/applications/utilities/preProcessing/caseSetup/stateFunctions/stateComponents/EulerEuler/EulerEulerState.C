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
    (c) 2016-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "EulerEulerState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "multiphaseThermo/multiphaseThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace stateFunctions
{

defineTypeNameAndDebug(EulerEulerState, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

EulerEulerState::EulerEulerState
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
        meshName,
        // Phasic turbulence
        true
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void EulerEulerState::initialise()
{
    for (const word& phaseName : phaseNames())
    {
        // fieldMaps
        stateDict_.subDict("fieldMaps").add
        (
            word("alpha." + phaseName),
            input().subDict("fieldMaps").lookup("alpha")
        );
        stateDict_.subDict("fieldMaps").add
        (
            word("U." + phaseName),
            input().subDict("fieldMaps").lookup("U")
        );
        if (input().subDict("fieldMaps").found("T"))
        {
            stateDict_.subDict("fieldMaps").add
            (
                word("T." + phaseName),
                input().subDict("fieldMaps").lookup("T")
            );
        }
    }
    // Remove old one from input to prevent it being merged later
    input().subDict("fieldMaps").remove("alpha");
    input().subDict("fieldMaps").remove("U");
    if (input().subDict("fieldMaps").found("T"))
    {
        input().subDict("fieldMaps").remove("T");
    }

    multiphaseMulticomponentState::initialise();

    // Check that we have a multiphase material model
    dictionary& matDict = system().subDict("materialProperties");
    const word matType(matDict.lookup("materialType"));
    if (matType != "multiphase")
    {
        FatalErrorInFunction
            << "The 'multiphaseEuler' state requires a multiphase material"
            << nl << exit(FatalError);
    }

    for (const word& phaseName : phaseNames())
    {
        // find grad(alpha) in defaults and write to state using proper phase name
        addAlphaTypeGrad(phaseName);
    }

    // Merge models from material defaults
    mergePhaseModels("phaseModel");

    // Merge binary models from defaults
    mergeBinaryPhaseModels("surfaceTensionModel", false);

    // Merge ordered binary models from defaults
    mergeBinaryPhaseModels("dragModel", true);
    mergeBinaryPhaseModels("dragModelSymmetric", false, false);
    mergeBinaryPhaseModels("liftModel", true, false);
    mergeBinaryPhaseModels("liftModelSymmetric", false, false);
    mergeBinaryPhaseModels("heatTransferModel", true, false);
    mergeBinaryPhaseModels("heatTransferModelSymmetric", false, false);
    mergeBinaryPhaseModels("massTransferModel", true, false);
    mergeBinaryPhaseModels("massTransferModelSymmetric", false, false);
    mergeBinaryPhaseModels("virtualMassModel", true, false);
    mergeBinaryPhaseModels("virtualMassModelSymmetric", false, false);
    mergeBinaryPhaseModels("wallLubrication", true, false);
    mergeBinaryPhaseModels("wallLubricationModelSymmetric", false, false);
    mergeBinaryPhaseModels("virtualMassModel", true, false);
    mergeBinaryPhaseModels("virtualMassModelSymmetric", false, false);
    mergeBinaryPhaseModels("aspectRatioModel", true, false);
    mergeBinaryPhaseModels("interfaceCompositionModel", true, false);
}


void EulerEulerState::finalise()
{
    multiphaseMulticomponentState::finalise();

    // Copy scheme to all phases
    dictionary& divSchemes(system().subDict("fvSchemes").subDict("divSchemes"));
    word oldScheme = "div(((rho*nuEff)*dev2(T(grad(U)))))";

    for (const word& phaseName : phaseNames())
    {
        divSchemes.add
        (
            word
            (
                "div((((alpha."+phaseName+"*thermo:rho."+phaseName
                +")*nuEff."+phaseName+")*dev2(T(grad(U."+phaseName+")))))"
            ),
            divSchemes.lookup(oldScheme),
            false
        );
    }
    divSchemes.remove(oldScheme);

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
