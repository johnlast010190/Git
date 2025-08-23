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

\*---------------------------------------------------------------------------*/

#include "db/error/error.H"

#include "submodels/WallInteractionModel/SpecularReflection/SpecularReflection.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SpecularReflection<CloudType>::SpecularReflection
(
    const dictionary& dict,
    CloudType& cloud
)
:
    WallInteractionModel<CloudType>(cloud)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SpecularReflection<CloudType>::~SpecularReflection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::SpecularReflection<CloudType>::correct
(
    typename CloudType::parcelType& p,
    const wallPolyPatch& wpp
)
{
    vector& U = p.U();

    vector nw = p.areaNormal();
    nw /= mag(nw);

    scalar U_dot_nw = U & nw;

    if (U_dot_nw > 0.0)
    {
        U -= 2.0*U_dot_nw*nw;
    }
}


// ************************************************************************* //
