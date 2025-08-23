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
    (c) 2019 OpenCFD Ltd.
    (c) 2024 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "combineMultipointFields/combineMultipointFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/Time/Time.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace functionObjects
    {
        defineTypeNameAndDebug(combineMultipointFields, 0);
        addToRunTimeSelectionTable(functionObject, combineMultipointFields, dictionary);
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
template<class Type>
void Foam::functionObjects::combineMultipointFields::calculateCombinedField
(
    Type& compositeField,
    const word& fieldName
)
{
    forAll(operatingPoints_, opPoints)
    {
        word regionName = operatingPoints_[opPoints].second();
        word geometryName = operatingPoints_[opPoints].first();
        const auto& geometryRegistry = time_.lookupObject<objectRegistry>
        (
           geometryName
        );
        const auto& regionRegistry = geometryRegistry.lookupObject<objectRegistry>
        (
            regionName
        );
        const Type& regionField =
            regionRegistry.lookupObject<Type>(fieldName);
        const auto& fi = geometryRegistry.lookupObject<volScalarField>("fi");
        forAll(fi, cI)
        {
            if (fi[cI]< 0)
            {
                compositeField[cI] = regionField[cI];
            }
        }
    }
    compositeField.write();
}

void Foam::functionObjects::combineMultipointFields::calculateAndWrite
(
    const word &fieldName,
    const word& newFieldName
)
{
    if (operatingPoints_.size())
    {
        word regionName = operatingPoints_[0].second();
        word geometryName = operatingPoints_[0].first();
        const auto& geometryRegistry = time_.lookupObject<objectRegistry>
        (
            geometryName
        );
        const auto& regionRegistry = geometryRegistry.lookupObject<objectRegistry>
        (
            regionName
        );
        if (regionRegistry.foundObject<volScalarField>(fieldName))
        {
            volScalarField compositeField
            (
                IOobject
                (
                    newFieldName,
                    time_.timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimless, 0)
            );
            calculateCombinedField(compositeField, fieldName);
        }
        else if (regionRegistry.foundObject<volVectorField>(fieldName))
        {
            volVectorField compositeField
            (
                IOobject
                (
                    newFieldName,
                    time_.timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector(dimless, vector::zero)
            );
            calculateCombinedField(compositeField, fieldName);
        }
        else if (regionRegistry.foundObject<volTensorField>(fieldName))
        {
            volTensorField compositeField
            (
                IOobject
                (
                    newFieldName,
                    time_.timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedTensor(dimless, tensor::zero)
            );
            calculateCombinedField(compositeField, fieldName);
        }
        else if (regionRegistry.foundObject<volSymmTensorField>(fieldName))
        {
            volSymmTensorField compositeField
            (
                IOobject
                (
                    newFieldName,
                    time_.timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedSymmTensor("0", dimless, symmTensor::zero)
            );
            calculateCombinedField(compositeField, fieldName);
        }
        else if (regionRegistry.foundObject<volSphericalTensorField>(fieldName))
        {
            volSphericalTensorField compositeField
            (
                IOobject
                (
                    newFieldName,
                    time_.timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedSphericalTensor("0", dimless, sphericalTensor::zero)
            );
            calculateCombinedField(compositeField, fieldName);
        }
        else
        {
            WarningInFunction
                << "Field " << fieldName
                << " could not be found in the object registry."
                << " No volume report compiled." << endl;
        }
    }
}

void Foam::functionObjects::combineMultipointFields::combineAndWriteGeometry()
{
    volScalarField totalFi
    (
        IOobject
        (
            totalFiName_,
            time_.timeName(),
            obr_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
         ),
         mesh_,
         dimensionedScalar(dimless, 0)
    );
    forAll(operatingPoints_, opPoints)
    {
        word geometryName = operatingPoints_[opPoints].first();
        const auto& geometryRegistry = time_.lookupObject<objectRegistry>
        (
            geometryName
        );
        const auto& fi = geometryRegistry.lookupObject<volScalarField>("fi");

        forAll(fi, cI)
        {
            if (fi[cI]< 0)
            {
                totalFi[cI] += fi[cI];
            }
        }
        const volScalarField::Boundary& fiBoundary = fi.boundaryField();
        volScalarField::Boundary& toalFiBoundary = totalFi.boundaryFieldRef();
        forAll(mesh_.boundary(), patchi)
        {
            forAll(fiBoundary[patchi], facei)
            {
                if (fiBoundary[patchi][facei]< 0)
                {
                    toalFiBoundary[patchi][facei] += fiBoundary[patchi][facei];
                }
            }
        }
    }
    totalFi.write();
}
// * * * * * * * * * * * * *      Constructors         * * * * * * * * * * * //

Foam::functionObjects::combineMultipointFields::combineMultipointFields
(
    const Foam::word &name,
    const Foam::Time &runTime,
    const Foam::dictionary& dict
)
:
fvMeshFunctionObject(name, runTime, dict),
fields_(0),
newFieldNames_(0),
totalFiName_(dict.lookupOrDefault<word>("totalFiName", "fi")),
operatingPoints_()
{
    read(dict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::combineMultipointFields::~combineMultipointFields()
= default;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool Foam::functionObjects::combineMultipointFields::execute()
{
    return true;
}

bool Foam::functionObjects::combineMultipointFields::write()
{
    Info<< type() << " " << name() << ":" << nl
        << "writing combined fields "<<endl;

    forAll(fields_, fieldI)
    {
        word fieldName = fields_[fieldI];
        word newFieldName = newFieldNames_[fieldI];
        calculateAndWrite(fieldName, newFieldName);
    }

    combineAndWriteGeometry();

    return true;
}

bool Foam::functionObjects::combineMultipointFields::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    Info<< type() << " " << name() << ":" << nl;

    {
        wordList fieldList = dict.lookup("fields");
        fields_.setSize(fieldList.size());
        fields_ = fieldList;
    }
    if (dict.found("newFieldNames"))
    {
        wordList fieldList = dict.lookup("newFieldNames");
        newFieldNames_.setSize(fieldList.size());
        newFieldNames_ = fieldList;
    }
    else
    {
        newFieldNames_ = fields_;
    }

    operatingPoints_ = dict.lookup<List<Tuple2<word, word>>>("operatingPoints");

    return true;
}
