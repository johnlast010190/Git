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
    (c) 2011-2023 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "submodels/Reacting/CompositionModel/CompositionModel/CompositionModel.H"
#include "materialModels/baseModels/materialModels.H"
#include "materialModels/materialTables/materialTables.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CompositionModel<CloudType>::CompositionModel(CloudType& owner)
:
    CloudSubModelBase<CloudType>(owner),
    thermo_(owner.thermo()),
    phaseProps_()
{}


template<class CloudType>
Foam::CompositionModel<CloudType>::CompositionModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type),
    thermo_(owner.thermo()),
    phaseProps_
    (
        this->coeffDict().lookup("phases"),
        thermo_.carrier().species(),
        thermo_.liquids().components(),
        thermo_.solids().components()
    )
{}


template<class CloudType>
Foam::CompositionModel<CloudType>::CompositionModel
(
    const CompositionModel<CloudType>& cm
)
:
    CloudSubModelBase<CloudType>(cm),
    thermo_(cm.thermo_),
    phaseProps_(cm.phaseProps_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CompositionModel<CloudType>::~CompositionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const Foam::SLGThermo& Foam::CompositionModel<CloudType>::thermo() const
{
    return thermo_;
}


template<class CloudType>
const Foam::basicSpecieMixture&
Foam::CompositionModel<CloudType>::carrier() const
{
    return thermo_.carrier();
}


template<class CloudType>
const Foam::liquidMixtureProperties&
Foam::CompositionModel<CloudType>::liquids() const
{
    return thermo_.liquids();
}


template<class CloudType>
const Foam::solidMixtureProperties&
Foam::CompositionModel<CloudType>::solids() const
{
    return thermo_.solids();
}


template<class CloudType>
const Foam::phasePropertiesList&
Foam::CompositionModel<CloudType>::phaseProps() const
{
    return phaseProps_;
}


template<class CloudType>
Foam::label Foam::CompositionModel<CloudType>::nPhase() const
{
    return phaseProps_.size();
}


template<class CloudType>
const Foam::wordList& Foam::CompositionModel<CloudType>::phaseTypes() const
{
    // if only 1 phase, return the constituent component names
    if (phaseProps_.size() == 1)
    {
        return phaseProps_[0].names();
    }
    else
    {
        return phaseProps_.phaseTypes();
    }
}


template<class CloudType>
const Foam::wordList& Foam::CompositionModel<CloudType>::stateLabels() const
{
    return phaseProps_.stateLabels();
}


template<class CloudType>
const Foam::wordList&
Foam::CompositionModel<CloudType>::componentNames(const label phasei) const
{
    return phaseProps_[phasei].names();
}


template<class CloudType>
Foam::label Foam::CompositionModel<CloudType>::carrierId
(
    const word& cmptName,
    const bool allowNotFound
) const
{
    label id = thermo_.carrierId(cmptName);

    if (id < 0 && !allowNotFound)
    {
        FatalErrorInFunction
            << "Unable to determine global id for requested component "
            << cmptName << ". Available components are " << nl
            << thermo_.carrier().species()
            << abort(FatalError);
    }

    return id;
}


template<class CloudType>
Foam::label Foam::CompositionModel<CloudType>::localId
(
    const label phasei,
    const word& cmptName,
    const bool allowNotFound
) const
{
    label id = phaseProps_[phasei].id(cmptName);

    if (id < 0 && !allowNotFound)
    {
        FatalErrorInFunction
            << "Unable to determine local id for component " << cmptName
            << abort(FatalError);
    }

    return id;
}


template<class CloudType>
Foam::label Foam::CompositionModel<CloudType>::localToCarrierId
(
    const label phasei,
    const label id,
    const bool allowNotFound
) const
{
    label cid = phaseProps_[phasei].carrierIds()[id];

    if (cid < 0 && !allowNotFound)
    {
        FatalErrorInFunction
            << "Unable to determine global carrier id for phase "
            << phasei << " with local id " << id
            << abort(FatalError);
    }

    return cid;
}


template<class CloudType>
const Foam::scalarField& Foam::CompositionModel<CloudType>::Y0
(
    const label phasei
) const
{
    return phaseProps_[phasei].Y();
}


template<class CloudType>
Foam::tmp<Foam::scalarField> Foam::CompositionModel<CloudType>::X
(
    const label phasei,
    const scalarField& Y
) const
{
    const phaseProperties& props = phaseProps_[phasei];
    scalarField X(Y.size());
    scalar WInv = 0.0;
    switch (props.phase())
    {
        case phaseProperties::GAS:
        {
            forAll(Y, i)
            {
                label cid = props.carrierIds()[i];
                X[i] = Y[i]/thermo_.carrier().W(cid);
                WInv += X[i];
            }
            break;
        }
        case phaseProperties::LIQUID:
        {
            forAll(Y, i)
            {
                X[i] = Y[i]/thermo_.liquids().properties()[i].W();
                WInv += X[i];
            }
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Only possible to convert gas and liquid mass fractions"
                << abort(FatalError);
        }
    }

    tmp<scalarField> tfld = X/(WInv + ROOTVSMALL);
    return tfld;
}


template<class CloudType>
Foam::scalar Foam::CompositionModel<CloudType>::H
(
    const label phasei,
    const scalarField& Y,
    const scalar p,
    const scalar T
) const
{
    const phaseProperties& props = phaseProps_[phasei];
    scalar HMixture = 0.0;
    switch (props.phase())
    {
        case phaseProperties::GAS:
        {
            forAll(Y, i)
            {
                label cid = props.carrierIds()[i];
                HMixture += Y[i]*const_cast<SLGThermo&>(thermo_).ha(cid, p, T);
            }
            break;
        }
        case phaseProperties::LIQUID:
        {
            forAll(Y, i)
            {
                HMixture += Y[i]*thermo_.liquids().properties()[i].ha(p, T);
            }
            break;
        }
        case phaseProperties::SOLID:
        {
            forAll(Y, i)
            {
                HMixture += Y[i]*thermo_.solids().properties()[i].ha(T);
            }
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown phase enumeration" << abort(FatalError);
        }
    }

    return HMixture;
}


template<class CloudType>
Foam::scalar Foam::CompositionModel<CloudType>::hs
(
    const label phasei,
    const scalarField& Y,
    const scalar p,
    const scalar T
) const
{
    const phaseProperties& props = phaseProps_[phasei];
    scalar hsMixture = 0.0;
    switch (props.phase())
    {
        case phaseProperties::GAS:
        {
            forAll(Y, i)
            {
                label cid = props.carrierIds()[i];
                const scalar hs =
                    const_cast<SLGThermo&>(thermo_).hs(cid, p, T);
                hsMixture += Y[i]*hs;
            }
            break;
        }
        case phaseProperties::LIQUID:
        {
            forAll(Y, i)
            {
                hsMixture +=
                    Y[i]*(thermo_.liquids().properties()[i].hs(p, T));
            }
            break;
        }
        case phaseProperties::SOLID:
        {
            forAll(Y, i)
            {
                hsMixture += Y[i]*thermo_.solids().properties()[i].hs(T);
            }
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown phase enumeration"
                << abort(FatalError);
        }
    }

    return hsMixture;
}


template<class CloudType>
Foam::scalar Foam::CompositionModel<CloudType>::Hc
(
    const label phasei,
    const scalarField& Y,
    const scalar p,
    const scalar T
) const
{
    const phaseProperties& props = phaseProps_[phasei];
    scalar HcMixture = 0.0;
    switch (props.phase())
    {
        case phaseProperties::GAS:
        {
            forAll(Y, i)
            {
                label cid = props.carrierIds()[i];
                HcMixture += Y[i]*const_cast<SLGThermo&>(thermo_).hf(cid);
            }
            break;
        }
        case phaseProperties::LIQUID:
        {
            forAll(Y, i)
            {
                HcMixture += Y[i]*thermo_.liquids().properties()[i].hf();
            }
            break;
        }
        case phaseProperties::SOLID:
        {
            forAll(Y, i)
            {
                HcMixture += Y[i]*thermo_.solids().properties()[i].hf();
            }
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown phase enumeration"
                << abort(FatalError);
        }
    }

    return HcMixture;
}


template<class CloudType>
Foam::scalar Foam::CompositionModel<CloudType>::Cp
(
    const label phasei,
    const scalarField& Y,
    const scalar p,
    const scalar T
) const
{
    const phaseProperties& props = phaseProps_[phasei];
    scalar CpMixture = 0.0;
    switch (props.phase())
    {
        case phaseProperties::GAS:
        {
            forAll(Y, i)
            {
                label cid = props.carrierIds()[i];
                CpMixture += Y[i]*const_cast<SLGThermo&>(thermo_).Cp(cid, p, T);
            }
            break;
        }
        case phaseProperties::LIQUID:
        {
            forAll(Y, i)
            {
                CpMixture += Y[i]*thermo_.liquids().properties()[i].Cp(p, T);
            }
            break;
        }
        case phaseProperties::SOLID:
        {
            forAll(Y, i)
            {
                CpMixture += Y[i]*thermo_.solids().properties()[i].Cp();
            }
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown phase enumeration"
                << abort(FatalError);
        }
    }

    return CpMixture;
}


template<class CloudType>
Foam::scalar Foam::CompositionModel<CloudType>::L
(
    const label phasei,
    const scalarField& Y,
    const scalar p,
    const scalar T
) const
{
    const phaseProperties& props = phaseProps_[phasei];
    scalar LMixture = 0.0;
    switch (props.phase())
    {
        case phaseProperties::GAS:
        {
            if (debug)
            {
                WarningInFunction
                    << "No support for gaseous components" << endl;
            }
            break;
        }
        case phaseProperties::LIQUID:
        {
            forAll(Y, i)
            {
                LMixture += Y[i]*thermo_.liquids().properties()[i].hl(p, T);
            }
            break;
        }
        case phaseProperties::SOLID:
        {
            if (debug)
            {
                WarningInFunction
                    << "No support for solid components" << endl;
            }
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown phase enumeration"
                << abort(FatalError);
        }
    }

    return LMixture;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "submodels/Reacting/CompositionModel/CompositionModel/CompositionModelNew.C"

// ************************************************************************* //
