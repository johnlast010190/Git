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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "UnsortedMeshedSurface/UnsortedMeshedSurface.H"
#include "containers/Lists/ListOps/ListOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
Foam::autoPtr<Foam::UnsortedMeshedSurface<Face>>
Foam::UnsortedMeshedSurface<Face>::New(const fileName& name, const word& ext)
{
    if (debug)
    {
        InfoInFunction << "Constructing UnsortedMeshedSurface" << endl;
    }

    auto cstrIter = fileExtensionConstructorTable_().find(ext);

    if (cstrIter == fileExtensionConstructorTable_().end())
    {
        // No direct reader, delegate to parent if possible
        const wordHashSet& delegate = ParentType::readTypes();

        if (delegate.found(ext))
        {
            // Create indirectly
            autoPtr<UnsortedMeshedSurface<Face>> surf
            (
                new UnsortedMeshedSurface<Face>
            );
            surf().transfer(ParentType::New(name, ext)());

            return surf;
        }
        else
        {
            FatalErrorInFunction
                << "Unknown file extension " << ext << nl << nl
                << exit(FatalError);
        }
    }

    return autoPtr<UnsortedMeshedSurface<Face>>(cstrIter->second(name));
}


template<class Face>
Foam::autoPtr<Foam::UnsortedMeshedSurface<Face>>
Foam::UnsortedMeshedSurface<Face>::New(const fileName& name)
{
    word ext = name.ext();
    if (ext == "gz")
    {
        ext = name.lessExt().ext();
    }

    return New(name, ext);
}


// ************************************************************************* //
