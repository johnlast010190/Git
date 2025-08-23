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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2017 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "SemiImplicitSource.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

template<class Type>
const Foam::wordList Foam::fv::SemiImplicitSource<Type>::volumeModeTypeNames_
(
    IStringStream("(absolute specific)")()
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type> template<class RhoFieldType>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fv::SemiImplicitSource<Type>::sourceMatrix
(
    const dimensionSet& dSet,
    const VolField<Type>& psi,
    const RhoFieldType& rho,
    const label fieldi
) const
{
    DebugInFunction
        << "SemiImplicitSource<" << pTraits<Type>::typeName
        << ">::addSup for source " << name_ << endl;

    if (useActivationField_)
    {
        if (useSigmoid_)
        {
            this->createSigmoidActivationField();
        }
        else
        {
            this->createActivationField();
        }
    }
    const word timeName(mesh_.time().timeName());
    DimensionedField<Type, volMesh> Su
    (
        IOobject(name_ + psi.name() + "Su", timeName, mesh_),
        mesh_,
        dimensioned<Type>("zero", dSet/dimVolume, Zero),
        false
    );

    if (SuTimeSeriesList_(fieldi) != nullptr)
    {
        UIndirectList<Type>(Su, cells_) =
            SuTimeSeriesList_[fieldi].value(mesh_.time().timeOutputValue())
           /(VDash_ + SMALL);
    }

    DimensionedField<scalar, volMesh> Sp
    (
        IOobject(name_ + psi.name() + "Sp", timeName, mesh_),
        mesh_,
        dimensionedScalar(Su.dimensions()/psi.dimensions(), 0),
        false
    );

    if (SpTimeSeriesList_(fieldi) != nullptr)
    {
        UIndirectList<scalar>(Sp, cells_) =
            SpTimeSeriesList_[fieldi].value(mesh_.time().timeOutputValue())
           /(VDash_ + SMALL);
    }

    if (useActivationField_)
    {
        return Su*activationField_() + fvm::SuSp(Sp*activationField_(), psi);
    }

    return Su*rho + fvm::SuSp(Sp*rho, psi);
}


template<class Type>
void Foam::fv::SemiImplicitSource<Type>::createActivationField() const
{
    scalar Vtotal = 0;

    const vector& activationPoint(activation_->value(mesh_.time().value()));

    activationField_() = 0;
    forAll(activationField_(), celli)
    {
        const vector& localCoords(csys().localPosition(mesh_.C()[celli]));

        if
        (
            (
                absoluteDistance_[0]
              ? mag(localCoords.x()) < activationPoint.x()
              : (
                    ( (localCoords.x() < activationPoint.x()) && (localCoords.x() >=0) ) ||
                    ( (localCoords.x() > activationPoint.x()) && (localCoords.x() <0) )
                )
            )
         && (
                absoluteDistance_[1]
              ? mag(localCoords.y()) < activationPoint.y()
              : (
                    ( (localCoords.y() < activationPoint.y()) && (localCoords.y() >=0) ) ||
                    ( (localCoords.y() > activationPoint.y()) && (localCoords.y() <0) )
                )
            )
         && (
                absoluteDistance_[2]
              ? mag(localCoords.z()) < activationPoint.z()
              : (
                    ( (localCoords.z() < activationPoint.z()) && (localCoords.z() >=0) ) ||
                    ( (localCoords.z() > activationPoint.z()) && (localCoords.z() <0) )
                )
            )
        )
        {
            activationField_()[celli] = 1.0;
            Vtotal += mesh_.V()[celli];
        }
    }

    reduce(Vtotal, sumOp<scalar>());
    DebugInFunction<< "Cell volume of the activation field " << Vtotal << endl;

    VDash_ = Vtotal;
}


template<class Type>
void Foam::fv::SemiImplicitSource<Type>::createSigmoidActivationField() const
{
    scalar Vtotal = 0;

    const vector& activationPoint(activation_->value(mesh_.time().value()));

    activationField_() = 0;
    forAll(activationField_(), celli)
    {
        const vector& localCoords(csys().localPosition(mesh_.C()[celli]));

        vector localC;

        forAll(localC, i)
        {
            if (absoluteDistance_[i])
                localC[i] = activationPoint[i] - mag(localCoords[i]);
            else
            {
                if (activationPoint[i] >= 0)
                {
                    if (localCoords[i] >= 0)
                    {
                        localC[i] = activationPoint[i] - localCoords[i];
                    }
                    else
                    {
                        localC[i] = localCoords[i];
                    }
                }
                else
                {
                    if (localCoords[i] <= 0)
                    {
                        localC[i] = localCoords[i] - activationPoint[i];
                    }
                    else
                    {
                        localC[i] = -localCoords[i];
                    }
                }
            }
        }
        activationField_()[celli] =
            sigmoid(localC.x(), sigmoidFactor_)*
            sigmoid(localC.y(), sigmoidFactor_)*
            sigmoid(localC.z(), sigmoidFactor_);

        Vtotal += activationField_()[celli]*mesh_.V()[celli];
    }

    reduce(Vtotal, sumOp<scalar>());

    DebugInFunction
        << "Cell volume of the Sigmoid activation field " << Vtotal << endl;

    VDash_ = Vtotal;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::scalar Foam::fv::SemiImplicitSource<Type>::sigmoid
(
    const scalar value,
    const scalar factor
) const
{
    return 1/(1 + Foam::exp(-value*factor));
}


template<class Type>
bool Foam::fv::SemiImplicitSource<Type>::wordListFind
(
    const wordList& wl,
    const word& w
)
{
    forAll(wl, i)
    {
        if (wl[i] == w)
        {
            return true;
        }
    }

    return false;
}


template<class Type>
void Foam::fv::SemiImplicitSource<Type>::setTimeSeries
(
    const dictionary& dict,
    wordList& fieldNames
)
{
    const dictionary& SuRateDict = dict.subDict("injectionRateSu");
    const dictionary& SpRateDict = dict.subDict("injectionRateSp");

    // check Su and Sp lists for non-unique field names
    {
        // add all unique field names to filedNames
        // toc and iterators seem smart enough not
        // to use duplicate
        fieldNames.setSize(SuRateDict.toc().size());
        label i = 0;

        forAllConstIter(dictionary, SuRateDict, iter)
        {
            fieldNames[i] = iter().keyword();
            i++;
        }

        wordList SpFields;
        SpFields.setSize(SpRateDict.toc().size());
        label j = 0;

        // grab all unique field names in SpRateDcit
        forAllConstIter(dictionary, SpRateDict, iter)
        {
            SpFields[j] = iter().keyword();
            j++;
        }

        // merge the two to make one unique list
        forAll(SpFields,fI)
        {
            if (!wordListFind(fieldNames, SpFields[fI]))
            {
                fieldNames.append(SpFields[fI]);
            }
        }
    }

    SuTimeSeriesList_.setSize(fieldNames.size());
    SpTimeSeriesList_.setSize(fieldNames.size());

    Info<< "    Reading Su and Sp semi-implicit source time series data"
        << endl;

    forAll(fieldNames, i)
    {
        fileName lookupFileName;

        if (SuRateDict.found(fieldNames[i]))
        {
            SuTimeSeriesList_.set
            (
                i,
                Function1<Type>::New(fieldNames[i], SuRateDict)
            );

            Info<< "      Found Field " << fieldNames[i]
                << " Su time series data" << endl;
        }
        else
        {
            SuTimeSeriesList_.set(i, nullptr);
            Info<< "      Field " << fieldNames[i]
                << " does not have Su time series data, assumed as zero. "
                << endl;
        }

        if (SpRateDict.found(fieldNames[i]))
        {
            SpTimeSeriesList_.set
            (
                i,
                Function1<scalar>::New(fieldNames[i], SpRateDict)
            );

            Info<< "      Found Field " << fieldNames[i]
                << " Sp time series data" << endl;

        }
        else
        {
            SpTimeSeriesList_.set(i, nullptr);
            Info<< "      Field " << fieldNames[i] << " does not have"
                << " Sp time series data, assumed as zero. " << endl;
        }
    }

    applied_.setSize(fieldNames_.size(), false);
}


template<class Type>
typename Foam::fv::SemiImplicitSource<Type>::volumeModeType
Foam::fv::SemiImplicitSource<Type>::wordToVolumeModeType
(
    const word& vmtName
) const
{
    forAll(volumeModeTypeNames_, i)
    {
        if (vmtName == volumeModeTypeNames_[i])
        {
            return volumeModeType(i);
        }
    }

    FatalErrorInFunction
        << "Unknown volumeMode type " << vmtName
        << ". Valid volumeMode types are:" << nl << volumeModeTypeNames_
        << exit(FatalError);

    return volumeModeType(0);
}


template<class Type>
Foam::word Foam::fv::SemiImplicitSource<Type>::volumeModeTypeToWord
(
    const volumeModeType& vmtType
) const
{
    if (vmtType > volumeModeTypeNames_.size())
    {
        return "UNKNOWN";
    }
    else
    {
        return volumeModeTypeNames_[vmtType];
    }
}


template<class Type>
void Foam::fv::SemiImplicitSource<Type>::convertOldInterface
(
    const dictionary& dict,
    wordList& fieldNames
)
{
    const dictionary& injectionDict(dict.subDict("injectionRateSuSp"));
    fieldNames.setSize(injectionDict.toc().size());
    List<Tuple2<Type, scalar>> injectionRate;
    injectionRate.setSize(fieldNames.size());

    label i = 0;
    forAllConstIter(dictionary, injectionDict, iter)
    {
        fieldNames[i] = iter().keyword();
        injectionDict.lookup(iter().keyword()) >> injectionRate[i];
        i++;
    }

    dictionary SuRateDict("injectionRateSu");
    dictionary SpRateDict("injectionRateSp");

    forAll(fieldNames, i)
    {
        SuRateDict.add(fieldNames[i], injectionRate[i].first());
        SpRateDict.add(fieldNames[i], injectionRate[i].second());
    }

    dictionary& dictRef = const_cast<dictionary&>(dict);
    dictRef.add("injectionRateSu", SuRateDict);
    dictRef.add("injectionRateSp", SpRateDict);

    DebugInFunction << "Converting old-style dictionary: " << dict << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::SemiImplicitSource<Type>::SemiImplicitSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    volumeMode_(vmAbsolute),
    VDash_(1.0),
    SuTimeSeriesList_(),
    SpTimeSeriesList_(),
    useActivationField_(dict.found("activation") ? true : false),
    useSigmoid_
    (
        dict.lookupOrDefault<Switch>("useSigmoidActivationFunction", false)
    ),
    sigmoidFactor_(dict.lookupOrDefault<scalar>("sigmoidFactor", 1.0)),
    activationField_
    (
        useActivationField_
      ? new volScalarField
        (
            IOobject
            (
                "activationField_" +name_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, 0)
        )
      : nullptr
    ),
    activation_
    (
        dict.found("activation")
      ? Function1<vector>::New("activation", dict)
      : nullptr
    ),
    absoluteDistance_
    (
        dict.lookupOrDefault<FixedList<bool, 3>>
        (
            "absoluteDistanceInActivationField",
            {0, 0, 0}
        )
    ),
    coorSysPtr_(nullptr),
    coorFramePtr_(nullptr)
{
    read(dict);
    createCS(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fv::SemiImplicitSource<Type>::sourceFields
(
    wordList& fieldNames
)
{
    if
    (
        coeffs_.isDict("injectionRateSuSp")
     && !coeffs_.isDict("injectionRateSu")
    )
    {
        convertOldInterface(coeffs_, fieldNames);
    }
    setTimeSeries(coeffs_, fieldNames);
}


template<class Type>
void Foam::fv::SemiImplicitSource<Type>::addSup
(
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    eqn +=
        this->sourceMatrix
        (
            eqn.dimensions(),
            eqn.psi(),
            geometricOneField(),
            fieldi
        );
}


template<class Type>
void Foam::fv::SemiImplicitSource<Type>::addSup
(
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    if (debug)
    {
        Info<< "SemiImplicitSource<" << pTraits<Type>::typeName
            << ">::addSup for source " << name_ << endl;
    }
    if (densityScaling_)
    {
        eqn +=
            this->sourceMatrix
            (
                eqn.dimensions()/rho.dimensions(),
                eqn.psi(),
                rho,
                fieldi
            );
    }
    else
    {
        this->addSup(eqn, fieldi);
    }
}


template<class Type>
void Foam::fv::SemiImplicitSource<Type>::addSup
(
    fvBlockMatrix<Type>& eqn,
    const label fieldi
)
{
    eqn +=
        this->sourceMatrix
        (
            eqn.dimensionSets()[0],
            eqn.psi(),
            geometricOneField(),
            fieldi
        );
}


template<class Type>
void Foam::fv::SemiImplicitSource<Type>::addSup
(
    const volScalarField& rho,
    fvBlockMatrix<Type>& eqn,
    const label fieldi
)
{
    if (debug)
    {
        Info<< "SemiImplicitSource<" << pTraits<Type>::typeName
            << ">::addSup for source " << name_ << endl;
    }

    if (densityScaling_)
    {
        eqn +=
            this->sourceMatrix
            (
                eqn.dimensionSets()[0]/rho.dimensions(),
                eqn.psi(),
                rho,
                fieldi
            );
    }
    else
    {
        this->addSup(eqn, fieldi);
    }
}


template<class Type>
void Foam::fv::SemiImplicitSource<Type>::createCS(const dictionary& dict)
{
    if (dict.found("referenceFrame"))
    {
        coorFramePtr_ = coordinateFrame::lookupNew(mesh_, dict);
    }
    else if (dict.found("coordinateSystem"))
    {
        coorSysPtr_ =
            coordinateSystem::New(mesh_, dict, coordinateSystem::typeName_());
    }
    else
    {
        coorSysPtr_.reset(new coordSystem::cartesian());
    }
}


// ************************************************************************* //
