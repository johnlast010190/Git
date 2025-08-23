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
    (c) 2012 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::IOporosityModelList

Description
    List of porosity models with IO functionality

SourceFiles
    IOporosityModelList.C

\*---------------------------------------------------------------------------*/

#ifndef IOporosityModelList_H
#define IOporosityModelList_H

#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "cfdTools/general/porosityModel/porosityModel/porosityModelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class IOporosityModelList Declaration
\*---------------------------------------------------------------------------*/

class IOporosityModelList
:
    public IOdictionary,
    public porosityModelList
{
private:

    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh, bool) const;


public:

    // Constructors

        //- Construct from mesh
        IOporosityModelList
        (
            const fvMesh& mesh,
            bool readFromFvOptions = false
        );

        //- Construct from mesh and registry
        IOporosityModelList
        (
            const objectRegistry& obr,
            const fvMesh& mesh,
            bool readFromFvOptions = false
        );

        //- Disallow default bitwise copy construction
        IOporosityModelList(const IOporosityModelList&) = delete;


        //- Destructor
        virtual ~IOporosityModelList()
        {}


    // Member Functions

        //- Read dictionary
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const IOporosityModelList&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
