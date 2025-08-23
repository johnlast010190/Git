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
    And the precice adapter solver object is based on the preCICE-
    adapter for OpenFOAM.

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
    (c) 2023 Engys Ltd.

\*---------------------------------------------------------------------------*/

using namespace Foam;

template<class FieldType>
void preciceAdapter::Adapter::loadCheckpointFields
(
    const PtrList<FieldType>& storedFields
)
{
    for (const FieldType& fc : storedFields)
    {
        // Load the volume field
        FieldType& f = fc.db().template lookupObjectRef<FieldType>(fc.name());
        f.forceAssign(fc);

        int nOldTimes(f.nOldTimes());
        if (nOldTimes >= 1)
        {
            f.oldTime().forceAssign(fc.oldTime());
        }
        if (nOldTimes == 2)
        {
            f.oldTime().oldTime().forceAssign(fc.oldTime().oldTime());
        }
    }
}

template<class FieldType>
void preciceAdapter::Adapter::storeCheckpointFields
(
    PtrList<FieldType>& storedFields
)
{
    for (FieldType& f : storedFields)
    {
        f.forceAssign(f.db().template lookupObject<FieldType>(f.name()));
    }
}