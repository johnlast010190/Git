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
    (c) 2019 Engys Ltd.
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "fieldCalcs/surfaceFieldCalc/surfaceFieldCalc.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

inline const Foam::labelList&
Foam::functionObjects::fieldCalcs::surfaceFieldCalc::faceId() const
{
    return faceId_;
}


inline const Foam::labelList&
Foam::functionObjects::fieldCalcs::surfaceFieldCalc::facePatch() const
{
    return facePatchId_;
}


inline const Foam::boolList&
Foam::functionObjects::fieldCalcs::surfaceFieldCalc::faceFlip() const
{
    return faceFlip_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::functionObjects::fieldCalcs::surfaceFieldCalc::regionTypes&
Foam::functionObjects::fieldCalcs::surfaceFieldCalc::regionType() const
{
    return regionType_;
}


// ************************************************************************* //
