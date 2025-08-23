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

\*---------------------------------------------------------------------------*/

#include "straightLine3DDataSetProvider.H"

#include "dataStructs/objects/line3D/straightLine3DSourceData.H"

#include "engysLine3D.h"

namespace Foam::functionObjects::runTimeVis
{

StraightLine3DDataSetProvider::StraightLine3DDataSetProvider
(
        const std::string& name,
        const StraightLine3DSourceData& data
)
:
        Line3DDataSetProvider(name)
{
    engysLine3D* lineFilter = getFilter();
    lineFilter->FromStraightLine(data.point1.x(), data.point1.y(), data.point1.z(),
                                 data.point2.x(), data.point2.y(), data.point2.z(),
                                 data.samplingPoints);
}

} // End namespace Foam


// ************************************************************************* //
