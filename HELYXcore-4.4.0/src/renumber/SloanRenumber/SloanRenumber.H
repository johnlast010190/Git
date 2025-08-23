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
    (c) 2012-2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::SloanRenumber

Description
    Sloan renumbering algorithm

    E.g. http://www.apav.it/sito_ratio/file_pdf/ratio_4/capitolo_1.pdf

SourceFiles
    SloanRenumber.C

\*---------------------------------------------------------------------------*/

#ifndef SloanRenumber_H
#define SloanRenumber_H

#include "renumberMethod.H"
#include "primitives/bools/Switch/Switch.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class SloanRenumber Declaration
\*---------------------------------------------------------------------------*/

class SloanRenumber
:
    public renumberMethod
{
    // Private data

        const Switch reverse_;


public:

    //- Runtime type information
    TypeName("Sloan");


    // Constructors

        //- Construct given the renumber dictionary
        SloanRenumber(const dictionary& renumberDict);

        //- Disallow default bitwise assignment
        SloanRenumber(const SloanRenumber&) = delete;


    //- Destructor
    virtual ~SloanRenumber()
    {}


    // Member Functions

        //- Return the order in which cells need to be visited, i.e.
        //  from ordered back to original cell label.
        //  This is only defined for geometric renumberMethods.
        virtual labelList renumber(const pointField&) const
        {
            NotImplemented;
        }

        //- Return the order in which cells need to be visited, i.e.
        //  from ordered back to original cell label.
        //  Use the mesh connectivity (if needed)
        virtual labelList renumber
        (
            const polyMesh& mesh,
            const pointField& cc
        ) const;

        //- Return the order in which cells need to be visited, i.e.
        //  from ordered back to original cell label.
        //  The connectivity is equal to mesh.cellCells() except
        //  - the connections are across coupled patches
        virtual labelList renumber
        (
            const labelListList& cellCells,
            const pointField& cc
        ) const;


    // Member Operators

        //- Disallow default bitwise copy construction
        void operator=(const SloanRenumber&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
