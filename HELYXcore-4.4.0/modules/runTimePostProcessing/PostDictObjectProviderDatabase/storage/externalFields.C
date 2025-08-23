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
    (c) 2022-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "externalFields.H"

#include "postDict/postDictKeys.H"

#include "db/solutionInstanceRegistry/solutionInstanceRegistry.H"

#undef Log
#include "vtkLogger.h"

namespace Foam::functionObjects::runTimeVis
{

ExternalFields::ExternalFields() = default;

void ExternalFields::updateForTimeStep(label timeIndex)
{
    domainRanges_.clear();
    if (providers_.empty()) return;
    for (const foamField& field : colourFields_.getColorFields())
    {
        vtkLogF(INFO, "Updating range for %s", field.c_str());

        scalarMinMax range;
        for (ExternalItemDataSetProvider* provider : providers_)
        {
            provider->updateIfNecessary(timeIndex, -1);
            range += provider->getRangeForField(field);
        }
        reduce(range, minMaxOp<scalar>());

        domainRanges_.set(field, range);
    }
}

void ExternalFields::addColourField(const foamField& colourField)
{
    colourFields_.addField(colourField);
}

scalarMinMax ExternalFields::getDomainRangeForField(const foamField &field) const
{
    return domainRanges_(field.lessAssociation(), Foam::scalarMinMax());
}

void ExternalFields::insertProvider(ExternalItemDataSetProvider *p)
{
    providers_.push_back(p);
}

void ExternalFields::clearProviders()
{
    providers_.clear();
}
// ************************************************************************* //
} // End namespace Foam
