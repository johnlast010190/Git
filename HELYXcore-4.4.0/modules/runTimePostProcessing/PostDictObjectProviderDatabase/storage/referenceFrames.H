/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : Dev
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
    (c) 2023-2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::ReferenceFrames

Description
    contains all reference frames defined in the meshObjects dictionary

SourceFiles
    ReferenceFrames.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "dictionaries/dictionaries.H"
#include "dataStructs/referenceFrame/referenceFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                        Class ReferenceFrames Declaration
\*---------------------------------------------------------------------------*/
class ReferenceFrames
{
private:
    HashTable<ReferenceFrame, word> frames_;
    const ReferenceFrame global_;

public:
    //- Gets the reference frames from the dict
    explicit ReferenceFrames(const Dictionaries& dictionaries, const fvMesh& defaultRegion);

    [[nodiscard]] const ReferenceFrame& getReferenceFrame(const word& name) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
