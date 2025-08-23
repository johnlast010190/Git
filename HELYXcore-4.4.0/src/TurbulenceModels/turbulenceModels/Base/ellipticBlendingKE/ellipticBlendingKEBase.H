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
    (c) 2019 Engys Ltd.


Class
    Foam::RASModels::v2fBase

Group
    grpRASTurbulence

Description
    Abstract base-class for v2-f models to provide BCs access to the v2 and f
    fields.

See also
    Foam::RASModels::v2f

SourceFiles
    v2fBase.C

\*---------------------------------------------------------------------------*/

#ifndef ellipticBlendingKEBase_H
#define ellipticBlendingKEBase_H

#include "RAS/RASModel/RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                            Class v2fBase Declaration
\*---------------------------------------------------------------------------*/

class ellipticBlendingKEBase
{
public:

    //- Runtime type information
//    TypeName("ellipticBlendingKEBase");


    // Constructors

        ellipticBlendingKEBase()
        {}


    //- Destructor
    virtual ~ellipticBlendingKEBase()
    {}


    // Member Functions

        //- Return turbulence stress normal to streamlines
        virtual tmp<volScalarField> v2k() const = 0;

        //- Return the damping function
        virtual tmp<volScalarField> a() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
