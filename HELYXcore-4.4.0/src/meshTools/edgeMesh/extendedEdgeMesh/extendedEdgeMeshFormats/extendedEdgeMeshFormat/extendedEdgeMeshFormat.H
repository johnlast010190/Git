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
    (c) 2013-2017 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::fileFormats::extendedEdgeMeshFormat

Description
    Provide a means of reading/writing the single-file OpenFOAM
    extendedEdgeMesh format

SourceFiles
    extendedEdgeMeshFormat.C

\*---------------------------------------------------------------------------*/

#ifndef extendedEdgeMeshFormat_H
#define extendedEdgeMeshFormat_H

#include "edgeMesh/extendedEdgeMesh/extendedEdgeMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                    Class extendedEdgeMeshFormat Declaration
\*---------------------------------------------------------------------------*/

class extendedEdgeMeshFormat
:
    public extendedEdgeMesh
{
public:

    // Constructors

        //- Construct from file name
        extendedEdgeMeshFormat(const fileName&);

        //- Disallow default bitwise copy construction
        extendedEdgeMeshFormat(const extendedEdgeMeshFormat&) = delete;


    //- Destructor
    virtual ~extendedEdgeMeshFormat()
    {}


    // Member Functions

        //- Read from file
        virtual bool read(const fileName&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const extendedEdgeMeshFormat&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
