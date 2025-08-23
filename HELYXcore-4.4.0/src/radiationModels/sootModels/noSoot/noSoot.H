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
    (c) 2013-2019 OpenFOAM Foundation

Class
    Foam::radiationModels::sootModels::noSoot

Description
    noSoot

SourceFiles
    noSoot.C

\*---------------------------------------------------------------------------*/

#ifndef noSoot_H
#define noSoot_H


#include "sootModels/sootModel/sootModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace sootModels
{

/*---------------------------------------------------------------------------*\
                             Class noSoot Declaration
\*---------------------------------------------------------------------------*/
class noSoot
:
    public sootModel
{

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        noSoot(const dictionary& dict, const fvMesh& mesh, const word&);


    //- Destructor
    virtual ~noSoot();


    // Member Functions

        // Edit

            //- Main update/correction routine
            void correct();

        // Access

            //- Return Ysoot
            const volScalarField& soot() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sootModels
} // End namespace Foam
} // End namespace radiationModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
