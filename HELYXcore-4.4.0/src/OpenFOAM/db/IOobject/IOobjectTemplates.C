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
    (c) 2015-17 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/IOobject/IOobject.H"
#include "global/fileOperations/fileOperation/fileOperation.H"
#include "db/IOstreams/Pstreams/Pstream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::IOobject::typeHeaderOk
(
    const bool checkType,
    const bool search,
    const bool verbose
)
{
    bool ok = true;

    // Everyone check or just master
    bool masterOnly =
        typeGlobal<Type>()
     && (
            IOobject::fileModificationChecking == timeStampMaster
         || IOobject::fileModificationChecking == inotifyMaster
        );

    const fileOperation& fp = Foam::fileHandler();

    // Determine local status
    if (!masterOnly || Pstream::master())
    {
        fileName fName(typeFilePath<Type>(*this, search));

        ok = fp.readHeader(*this, fName, Type::typeName);
        if (ok && checkType && headerClassName_ != Type::typeName)
        {
            if (verbose)
            {
                WarningInFunction
                    << "unexpected class name " << headerClassName_
                    << " expected " << Type::typeName
                    << " when reading " << fName << endl;
            }

            ok = false;
        }
    }

    // If masterOnly make sure all processors know about it
    if (masterOnly)
    {
        Pstream::scatter(ok);
    }

    return ok;
}


template<class Type>
void Foam::IOobject::warnNoRereading() const
{
    if (readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningInFunction
            << Type::typeName << ' ' << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED but "
            << Type::typeName << " does not support automatic rereading."
            << endl;
    }
}


// ************************************************************************* //
