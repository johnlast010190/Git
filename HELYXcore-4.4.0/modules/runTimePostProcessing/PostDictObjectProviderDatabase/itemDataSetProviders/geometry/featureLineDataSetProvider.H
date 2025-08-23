/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.0.1
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
    (c) 2022 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::FeatureLineDataSetProvider

Description
    Provider for the feature line objects that are loaded from eMesh files

SourceFiles
    featureLineDataSetProvider.C

\*---------------------------------------------------------------------------*/

#pragma once

#include <utility>

#include "itemDataSetProviders/geometry/surfaceDataSetProvider.H"
#include "dataStructs/geometry/featureLineData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations

namespace Foam::functionObjects::runTimeVis
{


/*---------------------------------------------------------------------------*\
              Class FeatureLineDataSetProvider Declaration
\*---------------------------------------------------------------------------*/

class FeatureLineDataSetProvider : public SurfaceDataSetProvider
{
private:

    FeatureLineData data_;

protected:

    void initialiseSurface() override;

public:
    FeatureLineDataSetProvider(
        const std::string &name,
        FeatureLineData data,
        const SurfaceTransformationsData& transformationsData
    )
        : SurfaceDataSetProvider(name, transformationsData), data_(std::move(data))
    {};

    ~FeatureLineDataSetProvider() override = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
