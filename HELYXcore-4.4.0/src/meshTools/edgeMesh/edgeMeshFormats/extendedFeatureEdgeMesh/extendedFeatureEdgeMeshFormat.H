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
    Foam::fileFormats::extendedFeatureEdgeMeshFormat

Description
    Provide a means of reading extendedFeatureEdgeMesh as featureEdgeMesh

SourceFiles
    extendedFeatureEdgeMeshFormat.C

\*---------------------------------------------------------------------------*/

#ifndef extendedFeatureEdgeMeshFormat_H
#define extendedFeatureEdgeMeshFormat_H

#include "edgeMesh/edgeMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                        Class extendedFeatureEdgeMeshFormat Declaration
\*---------------------------------------------------------------------------*/

class extendedFeatureEdgeMeshFormat
:
    public edgeMesh
{
public:

    // Constructors

        //- Construct from file name
        extendedFeatureEdgeMeshFormat(const fileName&);

        //- Disallow default bitwise copy construction
        extendedFeatureEdgeMeshFormat
        (
            const extendedFeatureEdgeMeshFormat&
        ) = delete;


    // Selectors

        //- Read file and return surface
        static autoPtr<edgeMesh> New(const fileName& name)
        {
            return autoPtr<edgeMesh>
            (
                new extendedFeatureEdgeMeshFormat(name)
            );
        }


    //- Destructor
    virtual ~extendedFeatureEdgeMeshFormat()
    {}


    // Member Functions

        //- Read from file
        virtual bool read(const fileName&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const extendedFeatureEdgeMeshFormat&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
