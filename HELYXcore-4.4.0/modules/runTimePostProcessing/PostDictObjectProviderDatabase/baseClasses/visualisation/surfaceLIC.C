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
    (c) 2019 OpenCFD Ltd.
    (c) 2022-2023 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::SurfaceLIC

Description
    Struct that holds Surface LIC SurfaceLIC information about an object.

\*---------------------------------------------------------------------------*/

#include "surfaceLIC.H"

#include "postDict/postDictKeys.H"

namespace Foam::functionObjects::runTimeVis
{

void SurfaceLIC::readFromDict(const dictionary& dict, const SurfaceLIC& defaultValues)
{
    visible = dict.lookupOrDefault<bool>(surfaceLICKeys::SHOW_KEY, defaultValues.visible);
    if (visible)
    {
        stepSize = dict.lookupOrDefault<scalar>(surfaceLICKeys::STEP_SIZE_KEY, defaultValues.stepSize);
        numberOfSteps = dict.lookupOrDefault<label>(surfaceLICKeys::NUMBER_OF_STEPS_KEY, defaultValues.numberOfSteps);
        vectorField = dict.lookup<word>(surfaceLICKeys::VECTOR_FIELD_KEY);
    }
}

bool SurfaceLIC::isVisible() const
{
    return visible;
}

SurfaceLIC& SurfaceLIC::operator=(const SurfaceLIC& other)
{
    if(this != &other)
    {
        this->visible = other.visible;
        this->stepSize = other.stepSize;
        this->numberOfSteps = other.numberOfSteps;
        this->vectorField = other.vectorField;
    }
    return *this;
}

} // End namespace Foam
