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
    (c) 2022-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "externalItemDataSetProvider.H"

#include "itemDataSetProvider.H"

#include "vtkDataSet.h"
#include "vtkDataSetAttributes.h"
#include "vtkCellData.h"
#include "vtkPointData.h"


namespace Foam::functionObjects::runTimeVis
{

ExternalItemDataSetProvider::ExternalItemDataSetProvider(std::string name)
    : ItemDataSetProvider(name)
{}

static scalarMinMax calculateRangeForField(vtkDataSetAttributes *attrs, const foamField& field)
{
    vtkDataArray *scalars = attrs->GetArray(field.getFoamName().c_str());
    if (!scalars)
    {
        return scalarMinMax();
    }

    double ranges[2];
    scalars->GetRange(ranges, field.getComponentIndex());
    return scalarMinMax(ranges[0], ranges[1]);
}

scalarMinMax ExternalItemDataSetProvider::getRangeForField(const foamField& field) const
{
    scalarMinMax range;
    vtkDataSet *ds = vtkDataSet::SafeDownCast(output_);
    if (ds)
    {
        if (field.isCellAssociation())
        {
            range += calculateRangeForField(ds->GetCellData(), field);
        }
        else
        {
            range += calculateRangeForField(ds->GetPointData(), field);
        }
    }
    return range;
}

} // End namespace Foam


// ************************************************************************* //
