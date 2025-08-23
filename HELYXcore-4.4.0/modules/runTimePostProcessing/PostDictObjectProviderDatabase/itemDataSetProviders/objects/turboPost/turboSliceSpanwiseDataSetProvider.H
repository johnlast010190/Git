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
    (c) 2023 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::TurboSliceSpanwiseDataSetProvider

Description
    Provider for the turbo spanwise slice post-processing objects,
    which is a slice with a constant stream

SourceFiles
    turboSliceSpanwiseObjectData.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "itemDataSetProviders/itemDataSetProvider.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class engysParametricSlicer;

namespace Foam::functionObjects::runTimeVis
{
struct TurboSliceSpanwiseObjectData;

/*---------------------------------------------------------------------------*\
              Class TurboSliceSpanwiseDataSetProvider Declaration
\*---------------------------------------------------------------------------*/

class TurboSliceSpanwiseDataSetProvider : public ItemDataSetProvider
{
private:

protected:

    vtkSmartPointer<engysParametricSlicer> turboSlicer_;

    void update(scalar currentTime) override;

public:
    TurboSliceSpanwiseDataSetProvider(const std::string& name, const TurboSliceSpanwiseObjectData& dictData);
    ~TurboSliceSpanwiseDataSetProvider() override = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
