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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "energyInverseTemperature.H"
#include "thermo/energy/matAbsoluteEnthalpy/matAbsoluteEnthalpy.H"
#include "thermo/energy/matAbsoluteInternalEnergy/matAbsoluteInternalEnergy.H"
#include "thermo/energy/matSensibleEnthalpy/matSensibleEnthalpy.H"
#include "thermo/energy/matSensibleInternalEnergy/matSensibleInternalEnergy.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Private Static Data * * * * * * * * * * * * * //

const Foam::scalar Foam::energyInverseTemperature::tol_ = 1.0e-4;

const int Foam::energyInverseTemperature::maxIter_ = 100;


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(energyInverseTemperature, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        energyInverseTemperature,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::energyInverseTemperature::lookupEnergyFieldName()
{
    const word eName = materialsDict().lookup<word>("energy");
    if (eName == matSensibleInternalEnergy::typeName)
    {
        return matSensibleInternalEnergy::name();
    }
    else if (eName == matSensibleEnthalpy::typeName)
    {
        return matSensibleEnthalpy::name();
    }
    else if (eName == matAbsoluteInternalEnergy::typeName)
    {
        return matAbsoluteInternalEnergy::name();
    }
    else if (eName == matAbsoluteEnthalpy::typeName)
    {
        return matAbsoluteEnthalpy::name();
    }
    else
    {
        FatalErrorInFunction
            << "Couldn't guess energy type for "
            << "the temperature inverse function. " << nl
            << exit(FatalError);
    }

    return word::null;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::energyInverseTemperature::energyInverseTemperature
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
:
    materialModel(obr, dict, phaseName, specieName, name)
{
    sMod_.setSize(modelsEnumSize_);
    dep_.setSize(5);
    energyInverseTemperature::read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::energyInverseTemperature::~energyInverseTemperature()
{}


// * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

void Foam::energyInverseTemperature::updateTable(const word& modelName)
{
    const word mixtureName =
        materialTables_.tableName(phaseName_, specieName_);
    if (!materialTables_.sTable().found(mixtureName))
    {
        FatalErrorInFunction
            << "Temperature model requires: \"" << mixtureName
            << "\" model."
            << abort(FatalError);
    }

    const matScalarTable& models = materialTables_.sTable()[mixtureName];

    if (modelName == TModel::typeName)
    {
        if
        (
            models.found(heModel::typeName)
         && models.found(CpvModel::typeName)
        )
        {
            // Required dependencies
            const wordList dependencies
            ({
                heModel::typeName,
                CpvModel::typeName,
                limitModel::typeName
            });

            // Dependency indices
            const labelList depInds({he, Cpv, limit});

            // Are models required?
            const boolList req({true, true, false});

            // Create links to dependent models
            fill(modelName, TModel::typeName, dependencies, depInds, req, 0);
        }
        else
        {
            FatalErrorInFunction
                << modelName << " model requires: " << mixtureName
                << " " << heModel::typeName << " and " << CpvModel::typeName << nl
                << " available models are: " << models.toc() << nl
                << " in table: \"" << mixtureName << "\"" << nl
                << abort(FatalError);
        }
    }
    else if (modelName == ThaModel::typeName)
    {
        if (models.found(haModel::typeName) && models.found(CpModel::typeName))
        {
            // Required dependencies
            const wordList dependencies
            ({
                haModel::typeName,
                CpModel::typeName,
                limitModel::typeName
            });

            // Dependency indices
            const labelList depInds({ha, Cp, limit});

            // Are models required?
            const boolList req({true, true, false});

            // Create links to dependent models
            fill(modelName, ThaModel::typeName, dependencies, depInds, req, 1);
        }
        else
        {
            FatalErrorInFunction
                << modelName << " model requires: " << mixtureName
                << " " << haModel::typeName << " and " << CpModel::typeName << nl
                << " available models are: " << models.toc() << nl
                << " in table: \"" << mixtureName << "\"" << nl
                << abort(FatalError);
        }
    }
    else if (modelName == ThsModel::typeName)
    {
        if (models.found(hsModel::typeName) && models.found(CpModel::typeName))
        {
            // Required dependencies
            const wordList dependencies
            ({
                hsModel::typeName,
                CpModel::typeName,
                limitModel::typeName
            });

            // Dependency indices
            const labelList depInds({hs, Cp, limit});

            // Are models required?
            const boolList req({true, true, false});

            // Create links to dependent models
            fill(modelName, ThsModel::typeName, dependencies, depInds, req, 2);
        }
        else
        {
            FatalErrorInFunction
                << modelName << " model requires: " << mixtureName
                << " " << hsModel::typeName << " and " << CpModel::typeName << nl
                << " available models are: " << models.toc() << nl
                << " in table: \"" << mixtureName << "\"" << nl
                << abort(FatalError);
        }
    }
    else if (modelName == TeaModel::typeName)
    {
        if (models.found(eaModel::typeName) && models.found(CvModel::typeName))
        {
            // Required dependencies
            const wordList dependencies
            ({
                eaModel::typeName,
                CvModel::typeName,
                limitModel::typeName
            });

            // Dependency indices
            const labelList depInds({ea, Cv, limit});

            // Are models required?
            const boolList req({true, true, false});

            // Create links to dependent models
            fill(modelName, TeaModel::typeName, dependencies, depInds, req, 3);
        }
        else
        {
            FatalErrorInFunction
                << modelName << " model requires: " << mixtureName
                << " " << eaModel::typeName << " and " << CvModel::typeName << nl
                << " available models are: " << models.toc() << nl
                << " in table: \"" << mixtureName << "\"" << nl
                << abort(FatalError);
        }
    }
    else if (modelName == TesModel::typeName)
    {
        if (models.found(esModel::typeName) && models.found(CvModel::typeName))
        {
            // Required dependencies
            const wordList dependencies
            ({
                esModel::typeName,
                CvModel::typeName,
                limitModel::typeName
            });

            // Dependency indices
            const labelList depInds({es, Cv, limit});

            // Are models required?
            const boolList req({true, true, false});

            // Create links to dependent models
            fill(modelName, TesModel::typeName, dependencies, depInds, req, 4);
        }
        else
        {
            FatalErrorInFunction
                << modelName << " model requires: " << mixtureName
                << " " << esModel::typeName << " and " << CvModel::typeName << nl
                << " available models are: " << models.toc() << nl
                << " in table: \"" << mixtureName << "\"" << nl
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorInFunction
            << modelName << " model requires: " << mixtureName
            << " " << heModel::typeName << " and " << CpvModel::typeName << nl
            << " available models are: " << models.toc() << nl
            << " in table: \"" << mixtureName << "\"" << nl
            << abort(FatalError);
    }
}


Foam::baseModels<Foam::scalar>* Foam::energyInverseTemperature::castScalarModel
(
    const word& modelName
)
{
    castMaterial(modelName, T)
    castMaterial(modelName, Ths)
    castMaterial(modelName, Tha)
    castMaterial(modelName, Tea)
    castMaterial(modelName, Tes)

    return nullptr;
}


bool Foam::energyInverseTemperature::updateScalarField
(
    const word& fieldName,
    const volScalarField& volField
)
{
    if (fieldName == phasePropertyName("T", phaseName_))
    {
        TnonConst_ = const_cast<volScalarField*>(&volField);
        return true;
    }
    else if
    (
        fieldName == phasePropertyName(lookupEnergyFieldName(), phaseName_)
    )
    {
        he_ = &volField;
        return true;
    }
    return false;
}


Foam::scalar Foam::energyInverseTemperature::Cell
(
    const label celli,
    const label heFun,
    const label CpvFun
) const
{
    scalar& Test = TnonConst_->operator[](celli);
    if (T_->isConst())
    {
        return Tref_ + T_->constant().value();
    }
    else if ((Test + Tref_) < 0)
    {
        FatalErrorInFunction
            << "Negative initial temperature T0: "
            << (Test + Tref_) << " [K]."
            << abort(FatalError);
    }

    //- Restore temperature after use
    const scalar Trestore = Test;

    scalar Tnew = Test;
    scalar Ttol = (Test + Tref_)*tol_;
    int iter = 0;

    do
    {
        Test = Tnew;
        Tnew =
            Test
          - (sMod_[heFun][celli] - he_->operator[](celli))
           /sMod_[CpvFun][celli];

        if (iter++ > maxIter_)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << maxIter_
                << abort(FatalError);
        }

    } while (mag(Tnew - Test) > Ttol);

    if (sMod_(limit) != nullptr)
    {
        Tnew = sMod_[limit][celli];
    }

    //- Restore temperature in registry
    Test = Trestore;

    return Tnew;
}


Foam::tmp<Foam::scalarField> Foam::energyInverseTemperature::Patch
(
    const label patchi,
    const label heFun,
    const label CpvFun
) const
{
    scalarField& Test = TnonConst_->boundaryFieldRef()[patchi];
    if (TnonConst_->boundaryFieldRef()[patchi].empty())
    {
        return tmp<scalarField>(new scalarField());
    }
    const scalarField& he = he_->boundaryField()[patchi];
    if (T_->isConst())
    {
        return
            tmp<scalarField>
            (
                new scalarField(Test.size(), Tref_ + T_->constant().value())
            );
    }
    else if ((min(Test) + Tref_) < 0)
    {
        FatalErrorInFunction
            << "Negative initial temperature T0: "
            << (min(Test) + Tref_) << " [K]."
            << abort(FatalError);
    }

    // Store old value of temperature
    scalarField Told(Test);

    tmp<scalarField> tTnew(new scalarField(Test));
    scalarField& Tnew = tTnew.ref();
    const scalar Ttol((min(Test) + Tref_)*tol_);
    int iter = 0;

    do
    {
        Test = Tnew;
        Tnew =
            Test
          - (sMod_[heFun].boundaryField()[patchi] - he)
           /sMod_[CpvFun].boundaryField()[patchi];

        if (iter++ > maxIter_)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << maxIter_
                << abort(FatalError);
        }
    } while (max(mag(Tnew - Test)) > Ttol);

    if (sMod_(limit) != nullptr)
    {
        Tnew = sMod_[limit].boundaryField()[patchi];
    }

    // Restore temperature in registry
    Test = Told;
    return tTnew;
}


Foam::tmp<Foam::scalarField> Foam::energyInverseTemperature::Internal
(
    const label heFun,
    const label CpvFun
) const
{
    scalarField& Test = TnonConst_->primitiveFieldRef();
    if (!Test.size())
    {
        return tmp<scalarField>(new scalarField());
    }
    const scalarField& he = he_->primitiveField();
    if (T_->isConst())
    {
        return
            tmp<scalarField>
            (
                new scalarField(Test.size(), Tref_ + T_->constant().value())
            );
    }
    else if ((min(Test) + Tref_) < 0)
    {
        FatalErrorInFunction
            << "Negative initial temperature T0: "
            << (min(Test) + Tref_) << " [K]."
            << abort(FatalError);
    }

    // Store old value of temperature
    scalarField Told(Test);

    tmp<scalarField> tTnew(new scalarField(Test));
    scalarField& Tnew = tTnew.ref();
    const scalar Ttol((min(Test) + Tref_)*tol_);

    int iter = 0;
    do
    {
        Test = Tnew;
        Tnew =
            Test
          - (sMod_[heFun].primitiveField() - he)
           /sMod_[CpvFun].primitiveField();
        if (iter++ > maxIter_)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << maxIter_
                << abort(FatalError);
        }
    } while (max(mag(Tnew - Test)) > Ttol);

    if (sMod_(limit) != nullptr)
    {
        Tnew = sMod_[limit].primitiveField();
    }

    // Restore temperature in registry
    Test = Told;
    return tTnew;
}


Foam::scalar Foam::energyInverseTemperature::TCell(const label celli) const
{
    return Cell(celli, he, Cpv);
}


Foam::scalar Foam::energyInverseTemperature::TValue(scalar p, scalar T) const
{
    NotImplemented;
}


Foam::tmp<Foam::scalarField> Foam::energyInverseTemperature::TPatch
(
    const label patchi
) const
{
    return Patch(patchi, he, Cpv);
}


Foam::tmp<Foam::scalarField> Foam::energyInverseTemperature::TInternal() const
{
    return Internal(he, Cpv);
}


Foam::tmp<Foam::volScalarField>
Foam::energyInverseTemperature::TGeometric() const
{
    tmp<volScalarField> tTest(volScalarField::New("Test", *TnonConst_));
    volScalarField& Test = tTest.ref();
    Test.primitiveFieldRef() = Internal(he, Cpv);
    forAll(Test.boundaryField(), patchi)
    {
        Test.boundaryFieldRef()[patchi].forceAssign(Patch(patchi, he, Cpv));
    }
    return tTest;
}


Foam::scalar Foam::energyInverseTemperature::ThaCell(const label celli) const
{
    return Cell(celli, ha, Cp);
}


Foam::scalar Foam::energyInverseTemperature::ThaValue(scalar p, scalar T) const
{
    NotImplemented;
}


Foam::tmp<Foam::scalarField> Foam::energyInverseTemperature::ThaPatch
(
    const label patchi
) const
{
    return Patch(patchi, ha, Cp);
}


Foam::tmp<Foam::scalarField>
Foam::energyInverseTemperature::ThaInternal() const
{
    return Internal(ha, Cp);
}


Foam::tmp<Foam::volScalarField>
Foam::energyInverseTemperature::ThaGeometric() const
{
    tmp<volScalarField> tTest(volScalarField::New("Test", *TnonConst_));
    volScalarField& Test = tTest.ref();
    Test.primitiveFieldRef() = Internal(ha, Cp);
    forAll(Test.boundaryField(), patchi)
    {
        Test.boundaryFieldRef()[patchi].forceAssign(Patch(patchi, ha, Cp));
    }
    return tTest;
}


Foam::scalar Foam::energyInverseTemperature::ThsCell(const label celli) const
{
    return Cell(celli, hs, Cp);
}


Foam::scalar Foam::energyInverseTemperature::ThsValue(scalar p, scalar T) const
{
    NotImplemented;
}


Foam::tmp<Foam::scalarField> Foam::energyInverseTemperature::ThsPatch
(
    const label patchi
) const
{
    return Patch(patchi, hs, Cp);
}


Foam::tmp<Foam::scalarField>
Foam::energyInverseTemperature::ThsInternal() const
{
    return Internal(hs, Cp);
}


Foam::tmp<Foam::volScalarField>
Foam::energyInverseTemperature::ThsGeometric() const
{
    tmp<volScalarField> tTest(volScalarField::New("Test", *TnonConst_));
    volScalarField& Test = tTest.ref();
    Test.primitiveFieldRef() = Internal(hs, Cp);
    forAll(Test.boundaryField(), patchi)
    {
        Test.boundaryFieldRef()[patchi].forceAssign(Patch(patchi, hs, Cp));
    }
    return tTest;
}


Foam::scalar Foam::energyInverseTemperature::TeaCell(const label celli) const
{
    return Cell(celli, ea, Cp);
}


Foam::scalar Foam::energyInverseTemperature::TeaValue(scalar p, scalar T) const
{
    NotImplemented;
}


Foam::tmp<Foam::scalarField> Foam::energyInverseTemperature::TeaPatch
(
    const label patchi
) const
{
    return Patch(patchi, ea, Cp);
}


Foam::tmp<Foam::scalarField>
Foam::energyInverseTemperature::TeaInternal() const
{
    return Internal(ea, Cp);
}


Foam::tmp<Foam::volScalarField>
Foam::energyInverseTemperature::TeaGeometric() const
{
    tmp<volScalarField> tTest(volScalarField::New("Test", *TnonConst_));
    volScalarField& Test = tTest.ref();
    Test.primitiveFieldRef() = Internal(ea, Cp);
    forAll(Test.boundaryField(), patchi)
    {
        Test.boundaryFieldRef()[patchi].forceAssign(Patch(patchi, ea, Cp));
    }
    return tTest;
}


Foam::scalar Foam::energyInverseTemperature::TesCell(const label celli) const
{
    return Cell(celli, es, Cp);
}


Foam::scalar Foam::energyInverseTemperature::TesValue(scalar p, scalar T) const
{
    NotImplemented;
}


Foam::tmp<Foam::scalarField> Foam::energyInverseTemperature::TesPatch
(
    const label patchi
) const
{
    return Patch(patchi, es, Cp);
}


Foam::tmp<Foam::scalarField>
Foam::energyInverseTemperature::TesInternal() const
{
    return Internal(es, Cp);
}


Foam::tmp<Foam::volScalarField>
Foam::energyInverseTemperature::TesGeometric() const
{
    tmp<volScalarField> tTest(volScalarField::New("Test", *TnonConst_));
    volScalarField& Test = tTest.ref();
    Test.primitiveFieldRef() = Internal(es, Cp);
    forAll(Test.boundaryField(), patchi)
    {
        Test.boundaryFieldRef()[patchi].forceAssign(Patch(patchi, es, Cp));
    }
    return tTest;
}


bool Foam::energyInverseTemperature::read()
{
    TnonConst_ = lookupPtr<scalar>(phasePropertyName("T", phaseName_));
    initialise("T");  // Phase is automatically added
    Tref_ = T_->offset().value();
    he_ =
        lookupPtr<scalar>
        (
            phasePropertyName(lookupEnergyFieldName(), phaseName_)
        );

    return true;
}


// ************************************************************************* //
