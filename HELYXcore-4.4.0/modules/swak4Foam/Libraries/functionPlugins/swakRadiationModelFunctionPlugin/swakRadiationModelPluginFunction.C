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
    (c) 2024 Engys Ltd.

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "swakRadiationModelPluginFunction.H"
#include "FieldValueExpressionDriver.H"

#include "containers/HashTables/HashPtrTable/HashPtrTable.H"

#include "include/swakThermoTypes.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#ifdef FOAM_RADIATION_HAS_SOOTMODEL
#include "sootModels/sootModel/sootModel.H"
#endif
#include "absorptionEmissionModels/absorptionEmissionModel/absorptionEmissionModel.H"
#include "basicThermo/basicThermo.H"

namespace Foam
{

defineTypeNameAndDebug(swakRadiationModelPluginFunction, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

swakRadiationModelPluginFunction::swakRadiationModelPluginFunction
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

const radiationModel &swakRadiationModelPluginFunction::radiationInternal
(
    const fvMesh& reg
)
{
    static HashPtrTable<radiationModel> radiation_;

    if (reg.foundObject<radiationModel>("radiationProperties"))
    {
        if (debug)
        {
            Info<< "swakRadiationModelPluginFunction::radiationInternal: "
                << "already in memory" << endl;
        }
        // Somebody else already registered this
        return reg.lookupObject<radiationModel>("radiationProperties");
    }
    if (!radiation_.found(reg.name()))
    {
        if (debug)
        {
            Info<< "swakRadiationModelPluginFunction::radiationInternal: "
                << "not yet in memory for " << reg.name() << endl;
        }

        dimensionedScalar TRef = basicThermo::TRefIfFound(reg);

        // Create it ourself because nobody registered it
        radiation_.set
        (
            reg.name(),
            radiationModel::New
            (
                reg.lookupObject<volScalarField>("T"),
                TRef
            ).ptr()
        );

        Info<< "Created radiation model. Calculating to get values ..."
            << endl;
        radiation_[reg.name()]->correct();
    }

    return *(radiation_[reg.name()]);
}


const radiationModel& swakRadiationModelPluginFunction::radiation()
{
    return radiationInternal(mesh());
}


// * * * * * * * * * * * * * * * Concrete implementations * * * * * * * * * //

#define concreteRadiationFunction(funcName,resultType)                        \
class swakRadiationModelPluginFunction_##funcName                             \
:                                                                             \
    public swakRadiationModelPluginFunction                                   \
{                                                                             \
public:                                                                       \
    TypeName("swakRadiationModelPluginFunction_"#funcName);                   \
                                                                              \
    swakRadiationModelPluginFunction_##funcName                               \
    (                                                                         \
        const FieldValueExpressionDriver& parentDriver,                       \
        const word& name                                                      \
    )                                                                         \
    :                                                                         \
        swakRadiationModelPluginFunction                                      \
        (                                                                     \
            parentDriver,                                                     \
            name,                                                             \
            #resultType                                                       \
        )                                                                     \
    {}                                                                        \
                                                                              \
    void doEvaluation()                                                       \
    {                                                                         \
        result().setObjectResult                                              \
        (                                                                     \
            autoPtr<resultType>                                               \
            (                                                                 \
                new resultType(#funcName "_rad", radiation().funcName())      \
            )                                                                 \
        );                                                                    \
    }                                                                         \
};                                                                            \
defineTypeNameAndDebug(swakRadiationModelPluginFunction_##funcName, 0);       \
addNamedToRunTimeSelectionTable                                               \
(                                                                             \
    FieldValuePluginFunction,                                                 \
    swakRadiationModelPluginFunction_##funcName,                              \
    name,                                                                     \
    radiation_##funcName                                                      \
);

concreteRadiationFunction(Rp, volScalarField);

class swakRadiationModelPluginFunction_Ru
:
    public swakRadiationModelPluginFunction
{
public:

    TypeName("swakRadiationModelPluginFunction_Ru");

    swakRadiationModelPluginFunction_Ru
    (
        const FieldValueExpressionDriver& parentDriver,
        const word& name
    )
    :
        swakRadiationModelPluginFunction
        (
            parentDriver,
            name,
            "volScalarField"
        )
    {}


    void doEvaluation()
    {
        tmp<DimensionedField<scalar,volMesh>> truRad(radiation().Ru());
        const DimensionedField<scalar,volMesh>& ruRad = truRad();

        autoPtr<volScalarField> val
        (
            new volScalarField
            (
                IOobject
                (
                    "radiationRu",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar(ruRad.dimensions(), 0),
                "zeroGradient"
            )
        );
#ifdef FOAM_NO_DIMENSIONEDINTERNAL_IN_GEOMETRIC
        const_cast<scalarField&>(val->internalField().field())
#else
        val->dimensionedInternalField()
#endif
            = ruRad;

        result().setObjectResult(val);
    }
};


defineTypeNameAndDebug(swakRadiationModelPluginFunction_Ru, 0);
addNamedToRunTimeSelectionTable
(
    FieldValuePluginFunction,
    swakRadiationModelPluginFunction_Ru,
    name,
    radiation_Ru
);


class swakRadiationModelPluginFunction_radSource
:
    public swakRadiationModelPluginFunction
{
public:

    TypeName("swakRadiationModelPluginFunction_radSource");

    swakRadiationModelPluginFunction_radSource
    (
        const FieldValueExpressionDriver& parentDriver,
        const word& name
    )
    :
        swakRadiationModelPluginFunction
        (
            parentDriver,
            name,
            "volScalarField"
        )
    {}

    void doEvaluation()
    {
        const volScalarField T4
        (
            pow4
            (
                mesh().lookupObject<volScalarField>("T")
              + basicThermo::TRefIfFound(mesh())
            )
        );

        autoPtr<volScalarField> val
        (
            new volScalarField
            (
                IOobject
                (
                    "radSource",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                -T4*radiation().Rp(),
                "zeroGradient"
            )
        );

#ifdef FOAM_NO_DIMENSIONEDINTERNAL_IN_GEOMETRIC
        const_cast<scalarField&>(val().internalField().field())
#else
        val().internalField()
#endif
            += radiation().Ru();
        val().correctBoundaryConditions();

        result().setObjectResult(val);
    }
};


defineTypeNameAndDebug(swakRadiationModelPluginFunction_radSource, 0);
addNamedToRunTimeSelectionTable
(
    FieldValuePluginFunction,
    swakRadiationModelPluginFunction_radSource,
    name,
    radiation_radSource
);

#ifdef FOAM_HAS_SOOTMODEL
    // Values of the sootModel
class swakRadiationModelPluginFunction_soot
:
    public swakRadiationModelPluginFunction
{
public:
    TypeName("swakRadiationModelPluginFunction_soot");
    swakRadiationModelPluginFunction_soot
    (
        const FieldValueExpressionDriver& parentDriver,
        const word& name
    )
    :
        swakRadiationModelPluginFunction
        (
            parentDriver,
            name,
            "volScalarField"
        )
    {}

    void doEvaluation()
    {
        result().setObjectResult
        (
            autoPtr<volScalarField>
            (
                new volScalarField("soot_rad", radiation().soot().soot())
            )
        );
    }
};

defineTypeNameAndDebug(swakRadiationModelPluginFunction_soot, 0);
addNamedToRunTimeSelectionTable
(
    FieldValuePluginFunction,
    swakRadiationModelPluginFunction_soot,
    name,
    radiation_soot
);

#endif

#ifdef FOAM_RADIATION_HAS_PUBLIC_ABSORPTIONMODEL
// values of the absortption-Emission submodel
#define concreteAbsEmiRadiationFunction(funcName,resultType)                 \
class swakRadiationModelPluginFunction_AbsEmi_##funcName                     \
:                                                                            \
    public swakRadiationModelPluginFunction                                  \
{                                                                            \
    label band_;                                                             \
public:                                                                      \
    TypeName("swakRadiationModelPluginFunction_AbsEmi_"#funcName);           \
    swakRadiationModelPluginFunction_AbsEmi_##funcName                       \
    (                                                                        \
        const FieldValueExpressionDriver& parentDriver,                      \
        const word& name                                                     \
    )                                                                        \
    :                                                                        \
        swakRadiationModelPluginFunction                                     \
        (                                                                    \
            parentDriver,                                                    \
            name,                                                            \
            #resultType,                                                     \
            "band primitive label"                                           \
        ),                                                                   \
        band_(-1)                                                            \
    {}                                                                       \
                                                                             \
    void doEvaluation()                                                      \
    {                                                                        \
        result().setObjectResult                                             \
        (                                                                    \
            autoPtr<resultType>                                              \
            (                                                                \
                new resultType                                               \
                (                                                            \
                    #funcName "_absEmi_rad",                                 \
                    radiation().absorptionEmission().funcName(band_)         \
                )                                                            \
            )                                                                \
        );                                                                   \
    }                                                                        \
    void setArgument                                                         \
    (                                                                        \
        label index,                                                         \
        const label& band                                                    \
    )                                                                        \
    {                                                                        \
        assert(index == 0);                                                  \
        band_ = band;                                                        \
    }                                                                        \
};                                                                           \
defineTypeNameAndDebug(swakRadiationModelPluginFunction_AbsEmi_##funcName, 0);  \
addNamedToRunTimeSelectionTable                                                 \
(                                                                               \
    FieldValuePluginFunction,swakRadiationModelPluginFunction_AbsEmi_##funcName,\
    name,                                                                       \
    radiation_AbsEmi_##funcName                                                 \
);

concreteAbsEmiRadiationFunction(a, volScalarField);
concreteAbsEmiRadiationFunction(aCont, volScalarField);
concreteAbsEmiRadiationFunction(aDisp, volScalarField);
concreteAbsEmiRadiationFunction(e, volScalarField);
concreteAbsEmiRadiationFunction(eCont, volScalarField);
concreteAbsEmiRadiationFunction(eDisp, volScalarField);
concreteAbsEmiRadiationFunction(E, volScalarField);
concreteAbsEmiRadiationFunction(ECont, volScalarField);
concreteAbsEmiRadiationFunction(EDisp, volScalarField);
#endif

} // namespace

// ************************************************************************* //
