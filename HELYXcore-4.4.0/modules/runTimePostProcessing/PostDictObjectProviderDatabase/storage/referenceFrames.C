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

\*---------------------------------------------------------------------------*/

#include "referenceFrames.H"
#include "postDict/postDictKeys.H"

namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * *       Constructors        * * * * * * * * * * * //

ReferenceFrames::ReferenceFrames(const Dictionaries& dictionaries, const fvMesh& defaultRegion)
: global_(defaultRegion)
{
    const dictionary& meshObjects = dictionaries.getMeshObjectsDict();
    for (const word& name : meshObjects.toc())
    {
        if (ReferenceFrame::isReferenceFrameDict(meshObjects, name))
        {
            frames_.insert(name, ReferenceFrame(defaultRegion, name));
        }
    }
    const dictionary &referenceFrames = GET_OPTIONAL_DICTIONARY(
        dictionaries.getRtppDict(),
        postDictKeys::REFERENCE_FRAMES_SUBDICT_KEY
    );
    for (const word& name : referenceFrames.toc())
    {
        if (ReferenceFrame::isReferenceFrameDict(referenceFrames, name))
        {
            const dictionary &subDict = referenceFrames.subDict(name);
            frames_.insert(name, ReferenceFrame(defaultRegion, name, subDict));
        }
    }
}

const ReferenceFrame& ReferenceFrames::getReferenceFrame(const word& referenceFrameName) const {
    if (referenceFrameName == meshObjectsKeys::GLOBAL_COORDINATE_SYSTEM_KEY ||
        referenceFrameName == meshObjectsKeys::INVALID_COORDINATE_SYSTEM_KEY ||
        referenceFrameName.empty())
    {
        return global_;
    }
    return frames_[referenceFrameName];
}

// ************************************************************************* //
} // End namespace
