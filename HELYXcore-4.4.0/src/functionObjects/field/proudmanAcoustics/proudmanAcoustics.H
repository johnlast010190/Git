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
    Foam::functionObjects::proudmanAcoustics

Group
    grpFieldFunctionObjects

Description
    Calculates and writes the acoustic sources based on
    the Proudman's formula.
    The turbulence kinetic energy and dissipation fields are retrieved
    from the turbulence model.

Usage
    Example of function object specification:
    \verbatim
    proudmanAcoustics1
    {
        type        proudmanAcoustics;
        libs        ("libfieldFunctionObjects.so");
        rhoInf      1.0;
        alpha0      343;
        alphae      0.1;
        refAcPower  1e-12;
        writedB     true;
        ...
    }
    \endverbatim

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    proudmanAcoustics.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_proudmanAcoustics_H
#define functionObjects_proudmanAcoustics_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "fields/volFields/volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class proudmanAcoustics Declaration
\*---------------------------------------------------------------------------*/

class proudmanAcoustics
:
    public fvMeshFunctionObject,
    public writeFile
{

protected:

    // Protected data

        //- Freestream density
        scalar rhoInf_;

        //- Speed of sound
        scalar alpha0_;

        //- Formula's constant
        scalar alphae_;

        //- Write power in dB
        Switch writedB_;

        //- Reference acoustic power
        scalar refAcPower_;

    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Calculate the shear-stress
        void calcProudmansFormula
        (
            volScalarField& acPower,
            const volScalarField& kTurb,
            const volScalarField& eTurb
        );


public:

    //- Runtime type information
    TypeName("proudmanAcoustics");


    // Constructors

        //- Construct from Time and dictionary
        proudmanAcoustics
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        proudmanAcoustics(const proudmanAcoustics&) = delete;


    //- Destructor
    virtual ~proudmanAcoustics();


    // Member Functions

        //- Read the proudmanAcoustics data
        virtual bool read(const dictionary&);

        //- Calculate the wall shear-stress
        virtual bool execute();

        //- Write the wall shear-stress
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const proudmanAcoustics&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
