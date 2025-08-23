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
    (c) 2019-2024 Engys Ltd.

Class
    Foam::openPatch

Description
    Activate patch if meanField on the patch is greater than a given value.
    Used in applications with opening patches like doors, windows etc.

SourceFiles
    openPatch.C

\*---------------------------------------------------------------------------*/

#ifndef openPatch_H
#define openPatch_H

#include "patchSwitch/switches/GIBSwitch/GIBSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshGIBChangers
{

/*---------------------------------------------------------------------------*\
                          Class openPatch Declaration
\*---------------------------------------------------------------------------*/

class openPatch
:
    public GIBSwitch
{
    // Private Data

        //- If the mean value of the 'fieldName_' patch field is above this
        //  trigger limit, the GIB will get disabled.
        const scalar fieldLim_;

        //- Name of the field
        const word fieldName_;

        //- Side of the GIB (master/slave) in which the condition is going to be
        //  evaluated
        const word side_;

        //- Evaluated patch ID
        label patchID_;


    // Private Member Functions

        //- The patch will never get enables.
        //  It is allowed to switch only from disable to enable.
        virtual bool enableCondition()
        {
            return false;
        }

        virtual bool disableCondition();


public:

    //- Runtime type information
    TypeName("openPatch");


    // Constructors

        //- Construct from fvMesh and dictionary
        openPatch(const fvMesh&, const dictionary&);

        //- Disallow default bitwise copy construction
        openPatch(const openPatch&) = delete;


    //- Destructor
    virtual ~openPatch(){};


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const openPatch&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshGIBChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
