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
    (c) 2015-2019 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::Line3DSourceInfo

Description
    Container for the cutting surface infos

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "db/dictionary/dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
class engysLine3D;

namespace Foam::functionObjects::runTimeVis
{
struct Id;
class ItemDataSetProvider;
class ItemRequirements;

class Line3DSourceInfo
{
public:
    Line3DSourceInfo() = default;

    virtual void readDict(const dictionary& sliceTypeDict, const string& caseFolder) = 0;

    [[nodiscard]] virtual ItemDataSetProvider* createDataSetProvider(const std::string& name) const = 0;

    virtual void addExtraSources(std::vector<Id> &sources) const
    {};

    virtual void addItemRequirements(ItemRequirements& itemRequirements) const {};

    virtual void computeAndAddToHash(size_t& hash) const = 0;

    virtual bool operator==(const Line3DSourceInfo& other) const = 0;

    virtual ~Line3DSourceInfo() = default;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
