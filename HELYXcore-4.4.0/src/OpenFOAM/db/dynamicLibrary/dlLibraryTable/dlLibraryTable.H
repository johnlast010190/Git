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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::dlLibraryTable

Description
    A table of dynamically loaded libraries

SourceFiles
    dlLibraryTable.C

\*---------------------------------------------------------------------------*/

#ifndef dlLibraryTable_H
#define dlLibraryTable_H

#include "primitives/ints/label/label.H"
#include "containers/Lists/DynamicList/DynamicList.H"
#include "primitives/strings/lists/fileNameList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class dlLibraryTable Declaration
\*---------------------------------------------------------------------------*/

class dlLibraryTable
{
    // Private Member Data

        //- List of pointers to the loaded dynamic libraries
        DynamicList<void*> libPtrs_;

        //- List of names of the loaded dynamic libraries
        DynamicList<fileName> libNames_;


public:

    // Declare name of the class and its debug switch
    ClassName("dlLibraryTable");


    // Constructors

        //- Construct null
        dlLibraryTable();

        //- Construct from dictionary and name of 'libs' entry giving
        //  the libraries to load
        dlLibraryTable(const dictionary&, const word&);

        //- Disallow default bitwise copy construction
        dlLibraryTable(const dlLibraryTable&) = delete;


    //- Destructor
    ~dlLibraryTable();


    // Member Functions

        //- Open the named library, optionally with warnings if problems occur
        bool open(const fileName& libName, const bool verbose = true);

        //- Open the named libraries, optionally with warnings if problems occur
        bool open(const fileNameList& libNames, const bool verbose = true);

        //- Open all the libraries listed in the 'libsEntry' entry in the
        //  given dictionary if present
        bool open(const dictionary&, const word& libsEntry);

        //- Find the handle of the named library
        void* findLibrary(const fileName& libName);

        //- Close the named library, optionally with warnings if problems occur
        bool close(const fileName& name, const bool verbose = true);

        //- Return the names of the loaded libraries
        const DynamicList<fileName>& libNames() const
        {
            return libNames_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dlLibraryTable&) = delete;
};


//- Table of loaded dynamic libraries
extern dlLibraryTable libs;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
