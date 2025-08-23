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

#include "turbulenceModelState.H"
#include "db/dictionary/functionEntries/includeEtcEntry/includeEtcEntry.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace stateFunctions
{

defineTypeNameAndDebug(turbulenceModelState, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

primitiveEntry turbulenceModelState::divScheme
(
    const word& fieldName, const dictionary& dict
) const
{
    return
    {
        phasicTurbulence_
      ? word("div(alphaPhi."+IOobject::group(fieldName)+","+fieldName+")")
      : word("div(phi," + fieldName + ")"),
        dict.lookup("divSchemes")
    };
}

void turbulenceModelState::addTurbulenceSpecificSchemes
(
    const dictionary& turbDict, dictionary& fvSchemes
)
{
    const dictionary& specificDivSchemes
    (
        turbDict.subDict("settings").subDict("divSchemes")
    );
    fvSchemes.subDict("divSchemes").merge(specificDivSchemes, true);
    const dictionary& specificLaplacianSchemes
    (
        turbDict.subDict("settings").subDict("laplacianSchemes")
    );
    fvSchemes.subDict("laplacianSchemes").merge(specificLaplacianSchemes, true);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulenceModelState::turbulenceModelState
(
    const word& region,
    const dictionary& input,
    const dictionary& defaults,
    const stateFunction& master,
    const stateIndex& index,
    const word& meshName,
    const bool phasicTurbulence
)
:
    regionState(region, input, defaults, master, index, meshName),
    phasicTurbulence_(phasicTurbulence)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void turbulenceModelState::fvSolutionInject
(
    const dictionary& fieldDict,
    const dictionary& modelDict,
    dictionary& fvSolDict,
    const word& tfname
)
{

    fvSolDict.subDict("relaxationFactors").subDict("equations").add
    (
        primitiveEntry
        (
            tfname,
            modelDict.lookup("relaxationFactor")
        ),
        false
    );

    //store existing entry (could be user defined state)
    dictionary prev;
    if (fvSolDict.subDict("solvers").found(tfname))
    {
        prev = fvSolDict.subDict("solvers").subDict(tfname);
        fvSolDict.subDict("solvers").remove(tfname);
    }

    primitiveEntry tol("vTol", modelDict.lookup("vTol"));
    primitiveEntry reltol
    (
        "vRelTol",
        modelDict.lookup("vRelTol")
    );

    fvSolDict.subDict("solvers").add(tfname, dictionary());
    dictionary& solverDict
    (
        fvSolDict.subDict("solvers").subDict(tfname)
    );
    solverDict.add(tol);
    solverDict.add(reltol);

    //cfg file
    fileName solverCfg
    (
        "caseDicts/preProcessing/caseSetup/settings/matrixSolvers/"
        + word(fieldDict.lookup("solver")) + ".cfg"
    );

    functionEntries::includeEtcEntry::mergeFile(solverDict, solverCfg);

    //overlay potential content from custom state
    fvSolDict.subDict("solvers").subDict(tfname).merge(prev);

    if (matrix() == maCoupled)
    {
        if (modelDict.lookup("maxIter"))
        {
            primitiveEntry maxIter
            (
                "maxIter",
                modelDict.lookup("maxIter")
            );
            solverDict.add(maxIter);
        }
    }

    //remove variables
    fvSolDict.subDict("solvers").subDict(tfname).remove("vTol");
    fvSolDict.subDict("solvers").subDict(tfname).remove("vRelTol");
}

void turbulenceModelState::initialise()
{
    wordList phases(phasicTurbulence_ ? phaseNames() : wordList(1, word::null));
    turbulenceFields_.setSize(phases.size());
    forAll(phases, phasei)
    {
        const word& phase = phases[phasei];

        //set turbulentProperties dictionary
        const word turbPropName =
            IOobject::groupName("turbulenceProperties", phase);
        constant().add(turbPropName, dictionary());
        dictionary& turbProp(constant().subDict(turbPropName));

        word turbCat(turbulenceTypeNames_[turbulence()]);

        turbProp.add("simulationType", turbCat);

        //get turbulence model defaults
        const dictionary& turbDB(defaults().subDict("turbulenceProperties"));

        word turbModelType
        (
            input().lookupOrDefault<word>
            (
                IOobject::groupName("turbulenceModel", phase),
                input().lookupOrDefault<word>("turbulenceModel", "laminar")
            )
        );

        //workaround for laminar that doesnt require redoing everything
        //for singular model category
        if (turbModelType == "laminar")
        {
            turbProp.add("simulationType", turbModelType, true);
        }


        if (turbulence() != tuLam && turbModelType != "laminar")
        {
            //add to fieldMaps and inject turbulenceProperties dictionary

            word turbModelCategory
            (
                (
                    // We revert to non-usf state for phasic. This is a bit of
                    // a hack to accommodate the fact that  EulerEuler doesn't
                    // yet follow the normal USF rules w.r.t bounded schemes.
                    (usf() && !phasicTurbulence_)
                  ? "USF" : word(compressibilityTypeNames_[compressibility()])
                )
              + word(turbulenceTypeNames_[turbulence()])
            );

            dictionary turbModelDefaults
            (
                turbDB.subDict(turbModelCategory).subDict(turbModelType)
            );

            //update fieldMaps
            if (turbModelDefaults.found("fieldMaps"))
            {
                const dictionary& turbFieldMaps
                (
                    turbModelDefaults.subDict("fieldMaps")
                );

                dictionary phaseTurbFieldMaps;
                forAllConstIters(turbFieldMaps, iter)
                {
                    autoPtr<entry> newEntry = iter().clone();
                    newEntry().keyword() =
                        IOobject::groupName(newEntry().keyword(), phase);
                    phaseTurbFieldMaps.add(newEntry);
                }

                stateDict_.subDict("fieldMaps").merge(phaseTurbFieldMaps);

                turbulenceFields_[phasei].setSize(phaseTurbFieldMaps.size());
                label i = 0;
                forAllConstIters(turbFieldMaps, iter)
                {
                    turbulenceFields_[phasei][i++] = iter().keyword();
                }

                turbModelDefaults.remove("fieldMaps");
            }

            //inject default settings into the specific turb model dictionary
            turbProp.add
            (
                turbCat, turbDB.subDict("baseTurbTypeDict")
            );

            dictionary& turbPropType(turbProp.subDict(turbCat));

            turbPropType.merge(turbModelDefaults);

            //do not overwrite
            turbPropType.add
            (
                word(turbCat + word("Model")), turbModelType, false
            );

            // add delta coefs
            if (turbulence() == tuLES && turbPropType.found("delta"))
            {
                word deltaType(turbModelDefaults.lookup("delta"));
                turbPropType.add
                (
                    word(deltaType + "Coeffs"),
                    turbDB.subDict(deltaType + "Coeffs")
                );
            }
        }
    }

    stateFunction::initialise();
}

void turbulenceModelState::correct()
{
    wordList phases(phasicTurbulence_ ? phaseNames() : wordList(1, word::null));
    forAll(phases, phasei)
    {
        const word& phase = phases[phasei];

        //skip turbulent field settings injection for laminar
        word turbModelType
        (
            input().lookupOrDefault<word>
            (
                IOobject::groupName("turbulenceModel", phase),
                input().lookupOrDefault<word>("turbulenceModel", "laminar")
            )
        );

        if (turbModelType != "laminar")
        {

            //inject fvSchemes & fvSolution settings for turbulence fields
            const word turbPropName =
                IOobject::groupName("turbulenceProperties", phase);
            word turbCat(turbulenceTypeNames_[turbulence()]);
            dictionary& turbCatDict
            (
                constant().subDict(turbPropName).subDict(turbCat)
            );


            forAll(turbulenceFields_[phasei], i)
            {
                //
                word tfmember(turbulenceFields_[phasei][i]);
                word tfname(IOobject::groupName(tfmember, phase));

                //get field definition
                const dictionary& turbFieldDict
                (
                    fieldDefinitions().subDict(tfname)
                );

                //check if field is solved for (solver)
                if (turbFieldDict.found("solver"))
                {
                    //get settings
                    word timeTN(timeTypeNames_[time()]);
                    word matrixTN(matrixTypeNames_[matrix()]);

                    const dictionary& turbSettings
                    (
                        turbCatDict.subDict("settings").
                            subDict(matrixTN).subDict(timeTN)
                    );

                    //fvSchemes
                    dictionary& fvSchemes(system().subDict("fvSchemes"));

                    fvSchemes.subDict("divSchemes").add
                    (
                        divScheme(tfname, turbSettings), false
                    );

                    addTurbulenceSpecificSchemes(turbCatDict, fvSchemes);

                    word lapScheme(word::null);

                    if (phasicTurbulence_)
                    {
                        if (compressibility() == ctIncomp)
                        {
                            lapScheme = "laplacian(D"+tfmember+"Eff,"+tfname+")";
                        }
                        else if (compressibility() == ctComp)
                        {
                            lapScheme =
                                "laplacian(((alpha."+phase+"*thermo:rho."
                               +phase+")*D"+tfmember+"Eff),"+tfname+")";
                        }

                        // Transform scheme which is specified in kOmegaSST
                        word oldScheme =
                            "laplacian((rho*D"+tfmember+"Eff),"+tfmember+")";
                        if
                        (
                            turbCatDict.subDict("settings").subDict
                            (
                                "laplacianSchemes"
                            ).found(oldScheme)
                        )
                        {
                            fvSchemes.subDict("laplacianSchemes").add
                            (
                                primitiveEntry
                                (
                                    lapScheme,
                                    turbCatDict.subDict("settings").subDict
                                    (
                                        "laplacianSchemes"
                                    ).lookup(oldScheme)
                                ),
                                false
                            );
                            fvSchemes.subDict("laplacianSchemes").remove
                            (
                                oldScheme
                            );
                        }
                        else
                        {
                            fvSchemes.subDict("laplacianSchemes").add
                            (
                                primitiveEntry
                                (
                                    lapScheme,
                                    turbSettings.lookup("laplacianSchemes")
                                ),
                                false
                            );
                        }

                        if (compressibility() == ctIncomp)
                        {
                            lapScheme =
                                "laplacian((alpha."+phase+"*nuEff."+phase+"),U."
                               +phase+")";
                        }
                        else if (compressibility() == ctComp)
                        {
                            lapScheme =
                                "laplacian(((alpha."+phase+"*thermo:rho."+phase
                               +")*nuEff."+phase+"),U."+phase+")";
                        }
                        fvSchemes.subDict("laplacianSchemes").add
                        (
                            primitiveEntry
                            (
                                lapScheme,
                                turbSettings.lookup("laplacianSchemes")
                            ),
                            false
                        );

                        // Transform scheme for kOmegaSST curvature
                        oldScheme = "div(phi,symm(grad(U)))";
                        if
                        (
                            turbCatDict.subDict("settings").subDict
                            (
                                "divSchemes"
                            ).found(oldScheme)
                        )
                        {
                            word divScheme =
                                "div(phiv."+phase+",symm(grad(U."+phase+")))";
                            fvSchemes.subDict("divSchemes").add
                            (
                                primitiveEntry
                                (
                                    divScheme,
                                    turbCatDict.subDict("settings").subDict
                                    (
                                        "divSchemes"
                                    ).lookup(oldScheme)
                                ),
                                false
                            );
                            fvSchemes.subDict("divSchemes").remove(oldScheme);
                        }
                    }
                    else
                    {
                        if (compressibility() == ctIncomp)
                        {
                            lapScheme = "laplacian(D"+tfname+"Eff,"+tfname+")";
                        }
                        else if (compressibility() == ctComp)
                        {
                            lapScheme = "laplacian((rho*D"+tfname+"Eff),"+tfname+")";
                        }
                        fvSchemes.subDict("laplacianSchemes").add
                        (
                            primitiveEntry
                            (
                                lapScheme,
                                turbSettings.lookup("laplacianSchemes")
                            ),
                            false
                        );
                    }

                    if (turbSettings.found("convectionGradient"))
                    {
                        fvSchemes.subDict("gradSchemes").add
                        (
                            primitiveEntry
                            (
                                "turbulence",
                                turbSettings.lookup("convectionGradient")
                            ),
                            false
                        );
                    }

                    //field gradients
                    if (turbSettings.found("fieldGradient"))
                    {
                        fvSchemes.subDict("gradSchemes").add
                        (
                            primitiveEntry
                            (
                                word("grad("+tfname+")"),
                                turbSettings.lookup("fieldGradient")
                            ),
                            false
                        );
                    }
                    if (turbSettings.found("nRequired"))
                    {
                        if (turbSettings.lookup<bool>("nRequired"))
                        {
                            fvSchemes.subDict("wallDist").add
                            (
                                "nRequired", turbSettings.lookup("nRequired")
                            );
                        }
                    }
                    //fvSolution
                    dictionary& fvSolution(system().subDict("fvSolution"));

                    fvSolutionInject
                    (
                        turbFieldDict,
                        turbSettings,
                        fvSolution,
                        tfname
                    );

                    if (time() == ttTrans)
                    {
                        fvSolutionInject
                        (
                            turbFieldDict,
                            turbSettings.subDict("final"),
                            fvSolution,
                            tfname+"Final"
                        );
                    }
                }
            }

            // Add generic phase-inversion handling as default for phasic models
            if (phasicTurbulence_)
            {
                dictionary stabDict;
                stabDict.add("type", "phaseTurbulenceStabilisation");
                stabDict.add("phase", phase);
                stabDict.add("alphaInversion", 0.1);
                if (!system().found("fvOptions"))
                {
                    system().add("fvOptions", dictionary());
                }
                system().subDict("fvOptions").add
                (
                    IOobject::groupName("phaseTurbulenceStabilisation", phase),
                    stabDict,
                    true
                );
            }

            if (turbCatDict.subDict("settings").found("cache"))
            {
                const dictionary& cacheDict = turbCatDict.subDict("settings").
                    subDict("cache");
                dictionary& fvSolution(system().subDict("fvSolution"));
                if (fvSolution.found("cache"))
                {
                    fvSolution.merge(cacheDict, true);
                }
                else
                {
                    fvSolution.add("cache", dictionary());
                    fvSolution.subDict("cache").merge(cacheDict, true);
                }
            }


            //remove settings
            if (turbCatDict.found("settings"))
            {
                turbCatDict.remove("settings");
            }
        }
    }
    stateFunction::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// ************************************************************************* //
