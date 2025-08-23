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
    2012-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "swakRhoChemistryModelPluginFunction.H"
#include "FieldValueExpressionDriver.H"
#include "containers/HashTables/HashPtrTable/HashPtrTable.H"
#include "include/swakThermoTypes.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "include/swak.H"
#include "fluidMulticomponentThermo/fluidMulticomponentThermo.H"
#include "multiphaseThermo/multiphaseThermo.H"

namespace Foam
{

defineTypeNameAndDebug(swakRhoChemistryModelPluginFunction, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

swakRhoChemistryModelPluginFunction::swakRhoChemistryModelPluginFunction
(
    const FieldValueExpressionDriver& parentDriver,
    const word& name,
    const word& returnValueType,
    const string& spec
):
    FieldValuePluginFunction
    (
        parentDriver,
        name,
        returnValueType,
        spec
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const basicChemistryModel&
swakRhoChemistryModelPluginFunction::chemistryInternal
(
    const fvMesh& reg
)
{
    static HashPtrTable<basicChemistryModel> chemistry_;

    // TODO: SWAK requires update for now it is switched off so it will
    // not lead to incorrect behavior, the thermo lookups have to be checked
    // The work should probably be done only on customer request.
    // Expecially it isn't clear where to get phase name from.
    NotImplemented;

    if (reg.foundObject<basicChemistryModel>("chemistryProperties"))
    {
        if (debug)
        {
            Info<< "swakRhoChemistryModelPluginFunction::chemistryInternal: "
                << "already in memory" << endl;
        }

        // Somebody else already registered this
        return reg.lookupObject<basicChemistryModel>("chemistryProperties");
    }

    if (!chemistry_.found(reg.name()))
    {
        if (debug)
        {
            Info<< "swakRhoChemistryModelPluginFunction::chemistryInternal: "
                << "not yet in memory for " << reg.name() << endl;
        }

        // Create it ourself because nobody registered it
        // TODO: Phase name isn't considered here yet... Needs to be
        // fixed in the future in order to to make use of it.
        // No test case is available for this yet so the rest of
        // the implementation should be considered ones the need
        // arises.
        fluidMulticomponentThermo& thermo =
            refCast<fluidMulticomponentThermo>
            (
                multiphaseThermo::lookupOrCreate(reg, word::null)
            );
        chemistry_.set(reg.name(), basicChemistryModel::New(thermo).ptr());

        Info<< "Created chemistry model. Calculating to get values ..."
            << endl;

        chemistry_[reg.name()]->solve
        (
#ifdef FOAM_CHEMISTRYMODEL_SOLVE_NEEDS_TIME
            reg.time().value(),
#endif
            reg.time().deltaT().value()
        );
    }

    return *(chemistry_[reg.name()]);
}


void swakRhoChemistryModelPluginFunction::updateChemistry(const scalar dt)
{
    const_cast<basicChemistryModel&>
    (
        chemistry()
    ).solve
    (
#ifdef FOAM_CHEMISTRYMODEL_SOLVE_NEEDS_TIME
        mesh().time().value(),
#endif
        dt
    );
}


#ifdef FOAM_RR_ONLY_DIMENSIONED_FIELD
tmp<volScalarField> swakRhoChemistryModelPluginFunction::wrapDimField
(
    const DimensionedField<scalar, volMesh>& dimField
)
{
    tmp<volScalarField> result
    (
        new volScalarField
        (
            IOobject
            (
                dimField.name(),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar(dimField.name(), dimField.dimensions(), 0),
            "zeroGradient"
        )
    );
#ifdef FOAM_NO_DIMENSIONEDINTERNAL_IN_GEOMETRIC
    const_cast<scalarField&>(result->internalField().field())
#else
    result->internalField()
#endif
    = dimField;

    return result;
}
#endif


const Foam::basicChemistryModel&
swakRhoChemistryModelPluginFunction::chemistry()
{
    return chemistryInternal(mesh());
}


// * * * * * * * * * * * * * * * Concrete implementations * * * * * * * * * //

#define concreteChemistryFunction(funcName,resultType)                        \
class swakRhoChemistryModelPluginFunction_##funcName                          \
:                                                                             \
    public swakRhoChemistryModelPluginFunction                                \
{                                                                             \
public:                                                                       \
    TypeName("swakRhoChemistryModelPluginFunction_"#funcName);                \
    swakRhoChemistryModelPluginFunction_##funcName                            \
    (                                                                         \
        const FieldValueExpressionDriver& parentDriver,                       \
        const word& name                                                      \
    )                                                                         \
    :                                                                         \
        swakRhoChemistryModelPluginFunction                                   \
        (                                                                     \
            parentDriver,                                                     \
            name,                                                             \
            #resultType                                                       \
        )                                                                     \
    {}                                                                        \
    void doEvaluation()                                                       \
    {                                                                         \
        result().setObjectResult                                              \
        (                                                                     \
            autoPtr<resultType>(new resultType(chemistry().funcName()))       \
        );                                                                    \
    }                                                                         \
};                                                                            \
defineTypeNameAndDebug(swakRhoChemistryModelPluginFunction_##funcName,0);     \
addNamedToRunTimeSelectionTable                                               \
(                                                                             \
    FieldValuePluginFunction,                                                 \
    swakRhoChemistryModelPluginFunction_##funcName,                           \
    name,                                                                     \
    psiChem_##funcName                                                        \
);

concreteChemistryFunction(tc,volScalarField);
#ifdef  FOAM_CHEMISTRYMODEL_HAS_NO_SOURCE_TERM
concreteChemistryFunction(Qdot,volScalarField);
#else
concreteChemistryFunction(Sh,volScalarField);
concreteChemistryFunction(dQ,volScalarField);
#endif

class swakRhoChemistryModelPluginFunction_RR
:
    public swakRhoChemistryModelPluginFunction
{
    // Private Member Data

        //- The species name
        word speciesName_;


public:

    TypeName("swakRhoChemistryModelPluginFunction_RR");

    swakRhoChemistryModelPluginFunction_RR
    (
        const FieldValueExpressionDriver& parentDriver,
        const word& name
    )
    :
        swakRhoChemistryModelPluginFunction
        (
            parentDriver,
            name,
            "volScalarField",
            "speciesName primitive word"
        )
    {}

    void doEvaluation()
    {
        label specI =
            chemistry().thermo().composition().species()[speciesName_];

        result().setObjectResult
        (
            autoPtr<volScalarField>
            (
#ifdef FOAM_RR_ONLY_DIMENSIONED_FIELD
                wrapDimField(chemistry().RR()[specI]).ptr()
#else
                new volScalarField(chemistry().RR(specI))
#endif
            )
        );
    }

    void setArgument(label index, const word& s)
    {
        assert(index == 0);
        speciesName_ = s;
    }
};


defineTypeNameAndDebug(swakRhoChemistryModelPluginFunction_RR, 0);
addNamedToRunTimeSelectionTable
(
    FieldValuePluginFunction,
    swakRhoChemistryModelPluginFunction_RR,
    name,
    psiChem_RR
);


class swakRhoChemistryModelPluginFunction_RRError
:
    public swakRhoChemistryModelPluginFunction
{
public:

    TypeName("swakRhoChemistryModelPluginFunction_RRError");

    swakRhoChemistryModelPluginFunction_RRError
    (
        const FieldValueExpressionDriver& parentDriver,
        const word& name
    )
    :
        swakRhoChemistryModelPluginFunction
        (
            parentDriver,
            name,
            "volScalarField"
        )
    {}

    void doEvaluation()
    {
        autoPtr<volScalarField> pSum
        (
            new volScalarField
            (
#ifdef FOAM_RR_ONLY_DIMENSIONED_FIELD
                0*wrapDimField(chemistry().RR()[0])
#else
                0*chemistry().RR(0)
#endif
            )
        );

        volScalarField& summe = pSum();
        for
        (
            label specI = 0;
            specI < chemistry().thermo().composition().species().size();
            specI++
        )
        {
#ifdef FOAM_RR_ONLY_DIMENSIONED_FIELD
            summe += wrapDimField(chemistry().RR()[specI]);
#else
            summe += chemistry().RR(specI);
#endif
        }

        result().setObjectResult(pSum);
    }
};


defineTypeNameAndDebug(swakRhoChemistryModelPluginFunction_RRError, 0);
addNamedToRunTimeSelectionTable
(
    FieldValuePluginFunction,
    swakRhoChemistryModelPluginFunction_RRError,
    name,
    psiChem_RRError
);


class swakRhoChemistryModelPluginFunction_RRSumPositive
:
    public swakRhoChemistryModelPluginFunction
{
public:
    TypeName("swakRhoChemistryModelPluginFunction_RRSumPositive");

    swakRhoChemistryModelPluginFunction_RRSumPositive
    (
        const FieldValueExpressionDriver& parentDriver,
        const word& name
    )
    :
        swakRhoChemistryModelPluginFunction
        (
            parentDriver,
            name,
            "volScalarField"
        )
    {}


    void doEvaluation()
    {
        autoPtr<volScalarField> pSum
        (
            new volScalarField
            (
#ifdef FOAM_RR_ONLY_DIMENSIONED_FIELD
                0*wrapDimField(chemistry().RR()[0])
#else
                0*chemistry().RR(0)
#endif
            )
        );

        volScalarField& summe = pSum();
        for
        (
            label specI = 0;
            specI < chemistry().thermo().composition().species().size();
            specI++
        )
        {
#ifdef FOAM_RR_ONLY_DIMENSIONED_FIELD
            tmp<volScalarField> tRR(this->wrapDimField(chemistry().RR()[specI]));
            const volScalarField& RR = tRR();
#else
            const volScalarField& RR=chemistry().RR(specI);
#endif
            forAll(summe, cellI)
            {
                if (RR[cellI] > 0)
                {
                    summe[cellI] += RR[cellI];
                }
            }
        }
        summe.correctBoundaryConditions();

        result().setObjectResult(pSum);
    }
};


defineTypeNameAndDebug(swakRhoChemistryModelPluginFunction_RRSumPositive, 0);
addNamedToRunTimeSelectionTable
(
    FieldValuePluginFunction,
    swakRhoChemistryModelPluginFunction_RRSumPositive,
    name,
    psiChem_RRSumPositive
);

class swakRhoChemistryModelPluginFunction_updateChemistry
:
    public swakRhoChemistryModelPluginFunction
{
    scalar timeStep_;


public:

    TypeName("swakRhoChemistryModelPluginFunction_updateChemistry");

    swakRhoChemistryModelPluginFunction_updateChemistry
    (
        const FieldValueExpressionDriver& parentDriver,
        const word& name
    )
    :
        swakRhoChemistryModelPluginFunction
        (
            parentDriver,
            name,
            "volScalarField",
            "timestep primitive scalar"
        )
    {}

    void doEvaluation()
    {
        updateChemistry(timeStep_);

        result().setObjectResult
        (
            autoPtr<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "dummyField",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedScalar(dimless, 0),
                    "zeroGradient"
                )
            )
        );
    }

    void setArgument(label index, const scalar& s)
    {
        assert(index == 0);
        timeStep_ = s;
    }
};


defineTypeNameAndDebug(swakRhoChemistryModelPluginFunction_updateChemistry, 0);
addNamedToRunTimeSelectionTable
(
    FieldValuePluginFunction,
    swakRhoChemistryModelPluginFunction_updateChemistry,
    name,
    psiChem_updateChemistry
);


class swakRhoChemistryModelPluginFunction_deltaTChem
:
    public swakRhoChemistryModelPluginFunction
{
public:

    TypeName("swakRhoChemistryModelPluginFunction_deltaTChem");

    swakRhoChemistryModelPluginFunction_deltaTChem
    (
        const FieldValueExpressionDriver& parentDriver,
        const word& name
    )
    :
        swakRhoChemistryModelPluginFunction
        (
            parentDriver,
            name,
            "volScalarField"
        )
    {}

    void doEvaluation()
    {
#ifdef FOAM_DELTATCHEM_NOT_DIMENSIONED
        const scalarField& dtChem = chemistry().deltaTChem();
#else
        const DimensionedField<scalar,volMesh>& dtChem =
            chemistry().deltaTChem();
#endif

        autoPtr<volScalarField> val
        (
            new volScalarField
            (
                IOobject
                (
                    "deltaTChem",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
#ifdef FOAM_DELTATCHEM_NOT_DIMENSIONED
                dimensionedScalar(dimless, 0),
#else
                dimensionedScalar(dtChem.dimensions(), 0),
#endif
                "zeroGradient"
            )
        );
#ifdef FOAM_DELTATCHEM_NOT_DIMENSIONED
        val->internalField() = dtChem;
#else
#ifdef FOAM_NO_DIMENSIONEDINTERNAL_IN_GEOMETRIC
        const_cast<scalarField&>(val->internalField().field()) = dtChem;
#else
        val->dimensionedInternalField() = dtChem;
#endif
#endif

        result().setObjectResult(val);
    }
};


defineTypeNameAndDebug(swakRhoChemistryModelPluginFunction_deltaTChem, 0);
addNamedToRunTimeSelectionTable
(
    FieldValuePluginFunction,
    swakRhoChemistryModelPluginFunction_deltaTChem,
    name,
    psiChem_deltaTChem
);

} // namespace

// ************************************************************************* //
