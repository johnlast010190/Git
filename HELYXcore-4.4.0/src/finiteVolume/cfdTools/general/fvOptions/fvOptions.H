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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2019-2024 Engys Ltd.

Class
    Foam::fv::options

Description
    Finite-volume options

SourceFiles
    options.C

\*---------------------------------------------------------------------------*/

#ifndef options_H
#define options_H

#include "cfdTools/general/fvOptions/fvOptionList.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "memory/autoPtr/autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class options Declaration
\*---------------------------------------------------------------------------*/

class options
:
    public IOdictionary,
    public optionList
{
protected:

    // Protected data

        //- Dictionary defaults to merge with the dictionary that is read drom
        //  disk
        dictionary injectDict_;


    // Protected Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;

        //- Create IO object if dictionary is present
        IOobject createIOobject
        (
            const fvMesh& mesh,
            const objectRegistry& obr
        ) const;


public:

    // Declare name of the class and its debug switch
    ClassName("fvOptions");


    // Constructors

        //- Construct from components with list of field names
        options
        (
            const fvMesh& mesh,
            const dictionary& injectDict = dictionary(),
            bool init = true
        );

        //- Construct from components with list of field names
        options
        (
            const fvMesh& mesh,
            const objectRegistry& obr,
            const dictionary& injectDict = dictionary(),
            bool init = true
        );

        //- Disallow default bitwise copy construction
        options(const options&) = delete;

        //- Construct fvOptions and register to database if not present
        //  otherwise lookup and return
        static options& New
        (
            const fvMesh& mesh,
            const dictionary& injectDict = dictionary(),
            bool init = true
        );

        //- Construct fvOptions and register to database if not present
        //  otherwise lookup and return
        static options& New
        (
            const fvMesh& mesh,
            const objectRegistry& obr,
            const dictionary& injectDict = dictionary(),
            bool init = true
        );


    //- Destructor
    virtual ~options()
    {}


    // Member Functions

        //- Inherit read from optionList
        using optionList::read;

        //- Read dictionary
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const options&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
