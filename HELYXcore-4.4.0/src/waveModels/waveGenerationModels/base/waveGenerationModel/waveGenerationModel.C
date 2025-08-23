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
    (c) 2016-2017 OpenCFD Ltd.
    (c) 2015 IH-Cantabria
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "waveGenerationModels/base/waveGenerationModel/waveGenerationModel.H"
#include "global/constants/mathematical/mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(waveGenerationModel, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::waveGenerationModel::waveGenerationModel
(
    const dictionary& dict,
    const fvMesh& mesh,
    const polyPatch& patch,
    const bool readFields
)
:
    waveModel(dict, mesh, patch, false),
    waveHeight_(0),
    waveAngle_(0)
{
    if (readFields)
    {
        readDict(dict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveModels::waveGenerationModel::~waveGenerationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::waveModels::waveGenerationModel::readDict
(
    const dictionary& overrideDict
)
{
    if (waveModel::readDict(overrideDict))
    {
        activeAbsorption_ = lookup<bool>("activeAbsorption");
        waveHeight_ = lookup<scalar>("waveHeight");
        if (waveHeight_ < 0)
        {
            FatalIOErrorInFunction(*this)
                << "Wave height must be greater than zero.  Supplied"
                << " value waveHeight = " << waveHeight_
                << exit(FatalIOError);
        }

        waveAngle_ = lookup<scalar>("waveAngle");
        waveAngle_ *= mathematical::pi/180;

        return true;
    }

    return false;
}


void Foam::waveModels::waveGenerationModel::info(Ostream& os) const
{
    waveModel::info(os);

    os  << "    Wave height : " << waveHeight_ << nl
        << "    Wave angle  : " << 180/mathematical::pi*waveAngle_ << nl;
}


// ************************************************************************* //
