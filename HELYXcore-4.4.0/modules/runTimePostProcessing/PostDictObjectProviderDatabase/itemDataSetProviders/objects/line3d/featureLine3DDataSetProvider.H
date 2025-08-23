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
    (c) 2020-2024 Engys Ltd.

Description
    Provider for the line 3D objects made from a feature line object

SourceFiles
    FeatureLine3DDataSetProvider.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "itemDataSetProviders/objects/line3d/line3DDataSetProvider.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{
struct FeatureLine3DSourceData;

/*---------------------------------------------------------------------------*\
 *            Class FeatureLine3DDataSetProvider Declaration
\*---------------------------------------------------------------------------*/

class FeatureLine3DDataSetProvider : public Line3DDataSetProvider
{
public:
    explicit FeatureLine3DDataSetProvider(const std::string& name);
    ~FeatureLine3DDataSetProvider() override = default;

protected:
    void updateFilter(engysLine3D* filter) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
