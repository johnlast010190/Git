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
    (c) 2020-2024 Engys Ltd.

Description
	Define primitive solid object for identifying whether part of a surface
	patch is inside a solid object

\*---------------------------------------------------------------------------*/

#include "acousticAnalogy/solidObject.H"

namespace Foam
{
    defineTypeNameAndDebug(solidObject, 0);
    defineRunTimeSelectionTable(solidObject, word);
}

Foam::solidObject::solidObject()
{
	r_=1;
	name_="none";
	type_="cylinder";
}

Foam::autoPtr<Foam::solidObject>
Foam::solidObject::New(
                        const word& name,
                        const fvMesh& mesh,
                        const dictionary& dict
                   )
{
    word typeName(dict.lookup("type"));

    const auto ctor =
        ctorTableLookup
        (
            "solidObject type",
            wordConstructorTable_(),
            typeName
        );
    return autoPtr<solidObject>(ctor(name, mesh, dict));
}


Foam::solidObject::solidObject
(
	const word& name,
	const fvMesh& mesh,
	const dictionary& dict
)
:
name_(name)
{
	type_="cylinder";
	r_=1;
}
