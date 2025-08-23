/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) ICE Stroemungsfoschungs GmbH
    (c) 2024 Engys Ltd.

Class
    recalcPhiFunctionObject

Description
    Recalculate phi according to the current velocity and density

SourceFiles
    recalcPhiFunctionObject.C

Contributors/Copyright:
    2008-2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef recalcPhiFunctionObject_H
#define recalcPhiFunctionObject_H

#include "update/updateSimple/updateSimpleFunctionObject.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class recalcPhiFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class recalcPhiFunctionObject
:
    public updateSimpleFunctionObject
{
    // Private Data

        //- Name of the velocity field
        word UName_;

        //- Name of the phi field
        word phiName_;

        //- Name of the pressure field
        word pName_;

        //- Name of the density field (if needed)
        word rhoName_;

        //- Write the fields?
        bool writeFields_;

        //- Write the old fields?
        bool writeOldFields_;


protected:

    // Protected Member Functions

        //- Do the actual phi calculation
        void recalc();


public:

    //- Runtime type information
    TypeName("recalcPhi");


    // Constructors

        //- Construct from components
        recalcPhiFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        recalcPhiFunctionObject(const recalcPhiFunctionObject&) = delete;


    // Member Functions

        bool start();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const recalcPhiFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
