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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "submodels/Reacting/PhaseChangeModel/PhaseChangeModel/PhaseChangeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class CloudType>
const Foam::wordList Foam::PhaseChangeModel<CloudType>::
enthalpyTransferTypeNames
(
    IStringStream
    (
        "("
            "latentHeat "
            "enthalpyDifference"
        ")"
    )()
);


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class CloudType>
typename Foam::PhaseChangeModel<CloudType>::enthalpyTransferType
Foam::PhaseChangeModel<CloudType>::wordToEnthalpyTransfer(const word& etName)
const
{
    forAll(enthalpyTransferTypeNames, i)
    {
        if (etName == enthalpyTransferTypeNames[i])
        {
            return enthalpyTransferType(i);
        }
    }

    FatalErrorInFunction
        << "Unknown enthalpyType " << etName << ". Valid selections are:" << nl
        << enthalpyTransferTypeNames << exit(FatalError);

    return enthalpyTransferType(0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PhaseChangeModel<CloudType>::PhaseChangeModel
(
    CloudType& owner
)
:
    CloudSubModelBase<CloudType>(owner),
    enthalpyTransfer_(etLatentHeat),
    dMass_(0.0)
{}


template<class CloudType>
Foam::PhaseChangeModel<CloudType>::PhaseChangeModel
(
    const PhaseChangeModel<CloudType>& pcm
)
:
    CloudSubModelBase<CloudType>(pcm),
    enthalpyTransfer_(pcm.enthalpyTransfer_),
    dMass_(pcm.dMass_)
{}


template<class CloudType>
Foam::PhaseChangeModel<CloudType>::PhaseChangeModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type),
    enthalpyTransfer_
    (
        wordToEnthalpyTransfer(this->coeffDict().lookup("enthalpyTransfer"))
    ),
    dMass_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PhaseChangeModel<CloudType>::~PhaseChangeModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const typename Foam::PhaseChangeModel<CloudType>::enthalpyTransferType&
Foam::PhaseChangeModel<CloudType>::enthalpyTransfer() const
{
    return enthalpyTransfer_;
}


template<class CloudType>
Foam::scalar Foam::PhaseChangeModel<CloudType>::dh
(
    const label idc,
    const label idl,
    const scalar p,
    const scalar T
) const
{
    return 0.0;
}


template<class CloudType>
Foam::scalar Foam::PhaseChangeModel<CloudType>::TMax
(
    const scalar p,
    const scalarField& X
) const
{
    return GREAT;
}


template<class CloudType>
Foam::scalar Foam::PhaseChangeModel<CloudType>::Tvap(const scalarField& X) const
{
    return -GREAT;
}


template<class CloudType>
void Foam::PhaseChangeModel<CloudType>::addToPhaseChangeMass(const scalar dMass)
{
    dMass_ += dMass;
}


template<class CloudType>
void Foam::PhaseChangeModel<CloudType>::info(Ostream& os)
{
    const scalar mass0 = this->template getBaseProperty<scalar>("mass");
    const scalar massTotal = mass0 + returnReduce(dMass_, sumOp<scalar>());

    Info<< "    Mass transfer phase change      = " << massTotal << nl;

    if (this->writeTime())
    {
        this->setBaseProperty("mass", massTotal);
        dMass_ = 0.0;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "submodels/Reacting/PhaseChangeModel/PhaseChangeModel/PhaseChangeModelNew.C"

// ************************************************************************* //
