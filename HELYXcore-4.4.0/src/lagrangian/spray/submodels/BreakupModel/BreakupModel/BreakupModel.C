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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "submodels/BreakupModel/BreakupModel/BreakupModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BreakupModel<CloudType>::BreakupModel
(
    CloudType& owner
)
:
    CloudSubModelBase<CloudType>(owner),
    solveOscillationEq_(false),
    y0_(0.0),
    yDot0_(0.0),
    TABComega_(0.0),
    TABCmu_(0.0),
    TABtwoWeCrit_(0.0)
{}


template<class CloudType>
Foam::BreakupModel<CloudType>::BreakupModel
(
    const BreakupModel<CloudType>& bum
)
:
    CloudSubModelBase<CloudType>(bum),
    solveOscillationEq_(bum.solveOscillationEq_),
    y0_(bum.y0_),
    yDot0_(bum.yDot0_),
    TABComega_(bum.TABComega_),
    TABCmu_(bum.TABCmu_),
    TABtwoWeCrit_(bum.TABtwoWeCrit_)
{}


template<class CloudType>
Foam::BreakupModel<CloudType>::BreakupModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type,
    bool solveOscillationEq
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type),
    solveOscillationEq_(solveOscillationEq),
    y0_(this->coeffDict().template lookupOrDefault<scalar>("y0", 0.0)),
    yDot0_(this->coeffDict().template lookupOrDefault<scalar>("yDot0", 0.0)),
    TABComega_(8),
    TABCmu_(5),
    TABtwoWeCrit_(12)
{
    if (solveOscillationEq_ && dict.found("TABCoeffs"))
    {
        const dictionary coeffs(dict.subDict("TABCoeffs"));
        TABComega_ = coeffs.lookup<scalar>("Comega");
        TABCmu_ = coeffs.lookup<scalar>("Cmu");
        TABtwoWeCrit_ = 2*coeffs.lookup<scalar>("WeCrit");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BreakupModel<CloudType>::~BreakupModel()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "submodels/BreakupModel/BreakupModel/BreakupModelNew.C"

// ************************************************************************* //
