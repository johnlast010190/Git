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
    (c) 2015-2016 OpenFOAM Foundation
    (c) 2021-2024 Engys Ltd

Class
    Foam::pisoControl

Description
    Specialization of the pimpleControl class for PISO control.

\*---------------------------------------------------------------------------*/

#ifndef pisoControl_H
#define pisoControl_H

#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"

//- Declare that pisoControl will be used
#define PISO_CONTROL
#undef PIMPLE_CONTROL

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class pisoControl Declaration
\*---------------------------------------------------------------------------*/

class pisoControl
:
    public pimpleControl
{
public:

    // Static Data Members

        //- Run-time type information
        TypeName("PISO");


    // Constructors

        //- Construct from mesh and the name of control sub-dictionary
        pisoControl(fvMesh& mesh, const word& dictName="PISO");

        //- Construct from object registry mesh and the name of control sub-dict
        pisoControl
        (
            fvMesh& mesh, const objectRegistry& obr, const word& dictName="PISO"
        );

        //- Disallow default bitwise copy construction
        pisoControl(const pisoControl&) = delete;


    //- Destructor
    virtual ~pisoControl();


    // Member Functions

        //- Return true if in the final inner (PISO) iteration
        inline bool finalInnerIter() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pisoControl&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cfdTools/general/solutionControl/pisoControl/pisoControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
