/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) ICE Stroemungsfoschungs GmbH
    (c) 1991-2008 OpenCFD Ltd.
    (c) 2025 Engys Ltd.

Contributors/Copyright:
    2012-2013, 2015-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "swakThermophysicalPluginFunction.H"
#include "FieldValueExpressionDriver.H"

#include "containers/HashTables/HashPtrTable/HashPtrTable.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#ifdef FOAM_HAS_FLUIDTHERMO
#include "solidThermo/solidThermo.H"
#endif

namespace Foam {

#ifdef FOAM_HAS_FLUIDTHERMO
defineTemplateTypeNameAndDebug(swakThermophysicalPluginFunction<swakFluidThermoType>,0);
defineTemplateTypeNameAndDebug(swakThermophysicalPluginFunction<solidThermo>,0);
#endif

defineTemplateTypeNameAndDebug(swakThermophysicalPluginFunction<basicThermo>,0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
swakThermophysicalPluginFunction<ThermoType>::swakThermophysicalPluginFunction(
    const FieldValueExpressionDriver &parentDriver,
    const word &name,
    const word &returnValueType
):
    FieldValuePluginFunction(
        parentDriver,
        name,
        returnValueType,
        string("")
    )
{
}

  // this explicit instantiation is needed on Ubunut 14.04
  // the compiler of that crappy pseudo-Linux needs it (but nobody else)
#ifdef FOAM_HAS_FLUIDTHERMO
template class swakThermophysicalPluginFunction<swakFluidThermoType>;
template class swakThermophysicalPluginFunction<solidThermo>;
#endif
template class swakThermophysicalPluginFunction<basicThermo>;

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const ThermoType &swakThermophysicalPluginFunction<ThermoType>::thermoInternal(
    const fvMesh &reg
)
{
    static HashPtrTable<ThermoType> thermo_;

    if (reg.foundObject<ThermoType>(basicThermo::dictName)) {
        if (debug) {
            Info<< "swakThermophysicalPluginFunction::thermoInternal: "
                << "already in memory" << endl;
        }
        // Somebody else already registered this
        if (reg.foundObject<ThermoType>(basicThermo::dictName)) {
            return reg.lookupObject<ThermoType>(basicThermo::dictName);
        }
    }
    if (!thermo_.found(reg.name())) {
        if (debug) {
            Info<< "swakThermophysicalPluginFunction::thermoInternal: "
                << "not yet in memory for " << reg.name() << endl;
        }

        bool found=false;

        {
            // make sure it is gone before we create the object
            IOobject thermoHeader
            (
                basicThermo::dictName,
                reg.time().caseConstant(),
                reg,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );
            IOobject materialsHeader
            (
                basicThermo::matDictName,
                reg.time().caseSystem(),
                reg,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );
            const bool isMaterials(materialsHeader.typeHeaderOk<IOdictionary>(true));

            IOdictionary dict(isMaterials ? materialsHeader : thermoHeader);

            word thermoTypeName;
            if (dict.isDict("thermoType")) {
                const dictionary& thermoTypeDict(dict.subDict("thermoType"));

                // Construct the name of the thermo package from the components
                thermoTypeName =
                    word(thermoTypeDict.lookup("type")) + '<'
                    + word(thermoTypeDict.lookup("mixture")) + '<'
                    + word(thermoTypeDict.lookup("transport")) + '<'
                    + word(thermoTypeDict.lookup("thermo")) + '<'
                    + word(thermoTypeDict.lookup("equationOfState")) + '<'
                    + word(thermoTypeDict.lookup("specie")) + ">>,"
                    + word(thermoTypeDict.lookup("energy")) + ">>>";
            } else {
                if (dict.found("thermoType")) {
                    FatalErrorInFunction
                        << "thermoType specification isn't supported option." << nl
                        << exit(FatalError);
                }

                thermoTypeName = dict.lookup<word>("materialType");
                if (thermoTypeName == "fluid" && dict.found("species")) {
                    thermoTypeName = "reactingFluid";
                }
            }

            swakRhoThermoType::objectRegistryConstructorTable::iterator cstrIter =
                swakRhoThermoType::objectRegistryConstructorTable_().find(
                    thermoTypeName
                );
            if (cstrIter != swakRhoThermoType::objectRegistryConstructorTable_().end())
            {
                if (debug) {
                    Info<< thermoTypeName << " is a rhoThermo-type";
                }
                found=true;
            } else if (debug) {
                Info<< "No " << thermoTypeName << " in rhoThermo-types "
                    << endl;
            }
        }

        if (!found) {
#ifdef FOAM_BASIC_THERMO_HAS_NO_NEW
            FatalErrorIn("swakThermophysicalPluginFunction<ThermoType>::thermoInternal")
                << "This version of Foam has no basicThermo::New"
                    << endl
                    << exit(FatalError);
#else
            thermo_.set(
                reg.name(),
                ThermoType::New(reg).ptr()
            );
#endif
        }
        else
        {
            // Create it ourself because nobody registered it
            thermo_.set
            (
                reg.name(),
                swakRhoThermoType::New(reg).ptr()
            );
        }
    }

    return *(thermo_[reg.name()]);
}

#ifdef FOAM_HAS_FLUIDTHERMO
template<>
const solidThermo &swakThermophysicalPluginFunction<solidThermo>::thermoInternal(
    const fvMesh &reg
)
{
    static HashPtrTable<solidThermo> thermo_;

    if (reg.foundObject<solidThermo>(basicThermo::dictName))
    {
        if (debug)
        {
            Info<< "swakThermophysicalPluginFunction::thermoInternal: "
                << "already in memory" << endl;
        }
        // Somebody else already registered this
        return reg.lookupObject<solidThermo>(basicThermo::dictName);
    }
    if (!thermo_.found(reg.name()))
    {
        if (debug)
        {
            Info<< "swakThermophysicalPluginFunction::thermoInternal: "
                << "not yet in memory for " << reg.name() << endl;
        }

        {
            // make sure it is gone before we create the object
            IOobject thermoHeader
            (
                basicThermo::dictName,
                reg.time().caseConstant(),
                reg,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );
            IOobject materialsHeader
            (
                basicThermo::matDictName,
                reg.time().caseSystem(),
                reg,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );
            const bool isMaterials(materialsHeader.typeHeaderOk<IOdictionary>(true));

            IOdictionary dict(isMaterials ? materialsHeader : thermoHeader);

            word thermoTypeName;
            if (dict.isDict("thermoType"))
            {
                const dictionary& thermoTypeDict(dict.subDict("thermoType"));

                // Construct the name of the thermo package from the components
                thermoTypeName =
                    word(thermoTypeDict.lookup("type")) + '<'
                    + word(thermoTypeDict.lookup("mixture")) + '<'
                    + word(thermoTypeDict.lookup("transport")) + '<'
                    + word(thermoTypeDict.lookup("thermo")) + '<'
                    + word(thermoTypeDict.lookup("equationOfState")) + '<'
                    + word(thermoTypeDict.lookup("specie")) + ">>,"
                    + word(thermoTypeDict.lookup("energy")) + ">>>";
            }
            else
            {
                if (dict.found("thermoType"))
                {
                    FatalErrorInFunction
                        << "thermoType specification isn't supported option." << nl
                        << exit(FatalError);
                }

                thermoTypeName = dict.lookup<word>("materialType");
                if (thermoTypeName == "fluid" && dict.found("species"))
                {
                    thermoTypeName = "reactingFluid";
                }
            }

            solidThermo::objectRegistryConstructorTable::iterator cstrIter =
                solidThermo::objectRegistryConstructorTable_().find
                (
                    thermoTypeName
                );
            if (cstrIter != solidThermo::objectRegistryConstructorTable_().end())
            {
                if (debug)
                {
                    Info<< thermoTypeName << " is a solidThermo-type";
                }
            }
            else if (debug)
            {
                Info<< "No " << thermoTypeName << " in solidThermo-types "
                    << endl;
            }

        }

        // Create it ourself because nobody registered it
        thermo_.set(
            reg.name(),
            solidThermo::New(reg).ptr()
        );
    }

    return *(thermo_[reg.name()]);
}
#endif

template<class ThermoType>
const ThermoType &swakThermophysicalPluginFunction<ThermoType>::thermo()
{
    return thermoInternal(mesh());
}

// * * * * * * * * * * * * * * * Concrete implementations * * * * * * * * * //

#define concreteThermoFunction(funcName,resultType,tthermo)        \
class swakThermophysicalPluginFunction_ ## funcName                \
: public swakThermophysicalPluginFunction<tthermo>                 \
{                                                                  \
public:                                                            \
    TypeName("swakThermophysicalPluginFunction_" #funcName);       \
    swakThermophysicalPluginFunction_ ## funcName (                \
        const FieldValueExpressionDriver &parentDriver,            \
        const word &name                                           \
    ): swakThermophysicalPluginFunction<tthermo>(                  \
        parentDriver,                                              \
        name,                                                      \
        #resultType                                                \
    ) {}                                                           \
    void doEvaluation() {                                          \
        result().setObjectResult(                                  \
            autoPtr<resultType>(                                   \
                new resultType(                                    \
                    thermo().funcName()                            \
                )                                                  \
            )                                                      \
        );                                                         \
    }                                                              \
};                                                                 \
defineTypeNameAndDebug(swakThermophysicalPluginFunction_ ## funcName,0);  \
addNamedToRunTimeSelectionTable(FieldValuePluginFunction,swakThermophysicalPluginFunction_ ## funcName,name,thermo_ ## funcName);

concreteThermoFunction(p,volScalarField,basicThermo);
concreteThermoFunction(rho,volScalarField,basicThermo);
concreteThermoFunction(psi,volScalarField,swakFluidThermoType);
#ifdef FOAM_HAS_FLUIDTHERMO
concreteThermoFunction(he,volScalarField,basicThermo);
concreteThermoFunction(Kappa,volVectorField,solidThermo);
#else
concreteThermoFunction(h,volScalarField,basicThermo);
concreteThermoFunction(hs,volScalarField,basicThermo);
concreteThermoFunction(e,volScalarField,basicThermo);
#endif
concreteThermoFunction(hf,volScalarField,basicThermo);
concreteThermoFunction(T,volScalarField,basicThermo);
concreteThermoFunction(Cp,volScalarField,basicThermo);
concreteThermoFunction(Cv,volScalarField,basicThermo);
concreteThermoFunction(mu,volScalarField,swakFluidThermoType);
concreteThermoFunction(kappa,volScalarField,basicThermo);

#ifdef FOAM_HAS_FLUIDTHERMO
// concreteThermoFunction(gamma,volScalarField,basicThermo);
concreteThermoFunction(Cpv,volScalarField,basicThermo);
// concreteThermoFunction(kappa,volScalarField,basicThermo);
    // concreteThermoFunction(kappaEff,volScalarField,basicThermo);
    // concreteThermoFunction(alphaEff,volScalarField,basicThermo);
#endif

} // namespace

// ************************************************************************* //
