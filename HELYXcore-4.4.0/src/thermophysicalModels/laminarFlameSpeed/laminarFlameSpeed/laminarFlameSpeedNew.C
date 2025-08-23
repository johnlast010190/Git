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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "laminarFlameSpeed/laminarFlameSpeed.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::laminarFlameSpeed> Foam::laminarFlameSpeed::New
(
    const psiuMulticomponentThermo& ct
)
{
    // do not register the dictionary
    IOdictionary propDict
    (
        IOobject
        (
            "combustionProperties",
            ct.T().time().constant(),
            ct.T().db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    const word corrType(propDict.lookup("laminarFlameSpeedCorrelation"));

    Info<< "Selecting laminar flame speed correlation " << corrType << endl;

    const auto ctor =
        ctorTableLookup
        (
            "laminarFlameSpeed type",
            dictionaryConstructorTable_(),
            corrType
        );
    return autoPtr<laminarFlameSpeed>(ctor(propDict, ct));
}


// ************************************************************************* //
