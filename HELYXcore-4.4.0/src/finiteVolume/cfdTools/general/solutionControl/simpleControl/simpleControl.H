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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2015-2024 Engys Ltd.

Class
    Foam::simpleControl

Description
    SIMPLE control class to supply convergence information/checks for
    the SIMPLE loop.

\*---------------------------------------------------------------------------*/

#ifndef simpleControl_H
#define simpleControl_H

#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"

//- Declare that simpleControl will be used
#define SIMPLE_CONTROL

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class simpleControl Declaration
\*---------------------------------------------------------------------------*/

class simpleControl
:
    public solutionControl
{
protected:

    // Protected Data

        //- Initialised flag
        bool initialised_;


    // Protected Member Functions

        //- Return true if equation variables should be stored (storePrevIter)
        virtual bool storeVars() const;


public:

    // Static Data Members

        //- Run-time type information
        TypeName("SIMPLE");


    // Constructors

        //- Construct from mesh
        simpleControl(fvMesh& mesh);

        //- Construct from mesh and registry
        simpleControl(fvMesh& mesh, const objectRegistry& obr);

        //- Disallow default bitwise copy construction
        simpleControl(const simpleControl&) = delete;


    //- Destructor
    virtual ~simpleControl();


    // Member Functions

        //- Read controls from fvSolution dictionary
        void readControls();

        //- Loop loop
        virtual bool loop();

        //- Return true if all convergence checks are satisfied
        bool criteriaSatisfied();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const simpleControl&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
