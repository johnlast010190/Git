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
    (c) 2017-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fieldInit.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fieldInit>
Foam::fieldInit::New
(
    const fvMesh& mesh,
    const fvSolutionRegistry& localDb,
    const dictionary& fieldDict,
    const word& fN
)
{

    const dictionary& initDict(fieldDict.subDict("initialisation"));

    //find list of allowed field inits from dict defaults
    wordList allowedFieldInits
    (
        fieldDict.lookup("allowedFieldInitialisationMethods")
    );

    //read initialization method
    word fieldInitTypeName
    (
        initDict.lookup("type")
    );

    bool allowed = false;

    //check if method is allowed according to field defaults
    forAll(allowedFieldInits, i)
    {
        if (allowedFieldInits[i] == fieldInitTypeName)
        {
            allowed = true;
        }
    }

    //if not, set to type none
    if (!allowed)
    {
            FatalErrorInFunction
                << "Initialisation method "<<fieldInitTypeName<<" for " << fN
                << ", not allowed.\n "
                << " Allowed methods for this field are "<< allowedFieldInits
                << exit(FatalError);
    }

    const auto ctor =
        ctorTableLookup
        (
            "Initialization type",
            initMethodConstructorTable_(),
            fieldInitTypeName
        );
    return autoPtr<fieldInit>(ctor(mesh, localDb, fieldDict, fN));
}


// ************************************************************************* //
