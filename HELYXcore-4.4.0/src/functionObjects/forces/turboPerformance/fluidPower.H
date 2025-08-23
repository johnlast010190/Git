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
    (c) D. Boger, B. Lewis M. Auvinen, H. Nilsson
    (c) 2024 Engys Ltd.

Class
    Foam::fluidPower

Description
    Calculates the mechanical energy change across the system
    (i.e. the work output).

    Member function calcDEmHead() calculates and returns the change
    in mechanical energy across the system and the hydrodynamical head.

    Member function fluidPower::write() calls calcDEmHead() and writes the
    data into the file <time dir>/forces.dat

Author
    Mikko Auvinen

Modifications
    Ported to OpenFoam-1.6-ext using patch by David Bogar - Penn State, ARL

    Added capability to analyse both pumps and turbines using a flag in the
    controlDict by Bryan Lewis - Penn State

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    fluidPower.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fluidPower_H
#define functionObjects_fluidPower_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "forces/forces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class forces Declaration
\*---------------------------------------------------------------------------*/

class fluidPower
:
    public forces
{
public:

    // Used to be 'pressureViscous'
    // Tuple which contains the work output (.first()) and the
    //hydrodynamic head (.second())
    typedef Tuple2<scalar , scalar> dEmHead;  // dEmHead := difference in
                                              //mechanical energy & Head

protected:

    // Private data

        // Read from dictonary

            //- Fag to determin if pump or turbine (Bryan)
            Switch turbine_;

            //- Patches to integrate the inflow of mechanical energy over
            labelHashSet inletPatchSet_;

            //- Patches to integrate the outflow of mechanical energy over
            labelHashSet outletPatchSet_;

            //- Name of flux field
            word phiName_;

            //- Sector based rotor-stator domain
            //  Inlet sections
            const label inletSectors_;

            //  Outlet sections
            const label outletSectors_;

            // volumetric flux
            mutable scalar volflux_;


    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("fluidPower");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        fluidPower
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool baseClass = false
        );

        //- Disallow default bitwise copy construction
        fluidPower(const fluidPower&) = delete;


    //- Destructor
    virtual ~fluidPower();


    // Member Functions

        //- Read the mechanical energy data
        virtual bool read(const dictionary&);

        //- Calculate the work output and write
        virtual bool write();

        //- Execute, currently does nothing
        virtual bool execute();

        //- Calculate and return the work output and hydrodynamic head
        virtual dEmHead calcDEmHead(const Switch&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fluidPower&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
