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
    (c) 2012-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::word&
Foam::fv::interRegionHeatTransferModel::nbrRegionName() const
{
    return nbrRegionName_;
}


inline const Foam::meshToMesh&
Foam::fv::interRegionHeatTransferModel::meshInterp() const
{
    if (!meshInterpPtr_.valid())
    {
        FatalErrorInFunction
            << "Interpolation object not set"
            << abort(FatalError);
    }

    return meshInterpPtr_();
}


inline const Foam::volScalarField&
Foam::fv::interRegionHeatTransferModel::htc() const
{
    return  htc_;
}


inline const Foam::fv::interRegionHeatTransferModel&
Foam::fv::interRegionHeatTransferModel::nbrModel() const
{
    if (nbrModel_ == nullptr)
    {
        FatalErrorInFunction
            << "Neighbour model not set"
            << abort(FatalError);
    }

    return *nbrModel_;
}


inline Foam::fv::interRegionHeatTransferModel&
Foam::fv::interRegionHeatTransferModel::nbrModel()
{
    if (nbrModel_ == nullptr)
    {
        FatalErrorInFunction
            << "Neighbour model not set"
            << abort(FatalError);
    }

    return *nbrModel_;
}


// ************************************************************************* //
