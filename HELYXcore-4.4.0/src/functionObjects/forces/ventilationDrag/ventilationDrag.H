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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2010-2024 Engys Ltd.

Class
    Foam::functionObjects::ventilationDrag

Group
    grpForcesFunctionObjects

Description
    Calculates the ventilation drag for specified wallVelocity patches, all
    belonging to the same reference frame.


Usage
    Example of function object specification:
    \verbatim
    ventilationDrag1
    {
        type               ventilationDrag;
        patches            ( wheels_RL_COUPE-rear-wheels-left wheels_RL_COUPE-rear-rim-left ); // list of patches of type wallVelocity belonging to the same rotational reference frame (separated by whitespace)
        Uinf               ( 38.89 0.0 0.0 );        // freestream velocity (equal to vehicle velocity)
        rhoInf             1.293;                    // freestream density (air density)
        referenceArea      1.79296625;
        functionObjectLibs ("libforces.so" );        // name of the library that has to be loaded to access the function object
    }
    \endverbatim


SourceFiles
    ventilationDrag.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_ventilationDrag_H
#define functionObjects_ventilationDrag_H

#include "forceCoeffs/forceCoeffs.H"
#include "referenceFrames/coordinateFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    class forceCoeffs;

/*---------------------------------------------------------------------------*\
                           Class ventilationDrag Declaration
\*---------------------------------------------------------------------------*/

class ventilationDrag
:
    public forces
{
    // Private data

        //- Wheel velocity field
        word Uname_;

        //- Velocity magnitude
        scalar magUInf_;

        //- Area
        scalar Aref_;

        //- Current reference frame pointer
        const coordinateFrame* coorFramePtr_;

        //- ventilationDrag file ptr
        autoPtr<OFstream> ventilationDragFilePtr_;


public:

    //- Runtime type information
    TypeName("ventilationDrag");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        ventilationDrag
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        ventilationDrag(const ventilationDrag&) = delete;


    //- Destructor
    virtual ~ventilationDrag();


    // Member Functions

        //- Read the parameters
        virtual bool read(const dictionary&) override;

        //- Execute
        virtual bool execute() override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ventilationDrag&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //