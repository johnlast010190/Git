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
    (c) 2011-2022 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "solidMulticomponentThermo/solidMulticomponentThermo.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidMulticomponentThermo, 0);
    defineRunTimeSelectionTable(solidMulticomponentThermo, objectRegistry);
    defineRunTimeSelectionTable(solidMulticomponentThermo, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidMulticomponentThermo::implementation::implementation
(
    const objectRegistry& obr,
    const word& phaseName
)
{}


Foam::solidMulticomponentThermo::implementation::implementation
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName
)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidMulticomponentThermo>
Foam::solidMulticomponentThermo::New
(
    const objectRegistry& obr,
    const word& phaseName
)
{
    return basicThermo::New<solidMulticomponentThermo>(obr, phaseName);
}


Foam::autoPtr<Foam::solidMulticomponentThermo>
Foam::solidMulticomponentThermo::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName
)
{
    return basicThermo::New<solidMulticomponentThermo>(obr, dict, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidMulticomponentThermo::~solidMulticomponentThermo()
{}


Foam::solidMulticomponentThermo::implementation::~implementation()
{}


// ************************************************************************* //
