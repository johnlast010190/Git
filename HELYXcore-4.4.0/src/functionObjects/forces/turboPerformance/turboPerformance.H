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
    (c) 2025 Engys Ltd.

Class
    Foam::functionObjects::turboPerformance

Authors
    Mikko Auvinen

Modifications
    Ported to OpenFoam-1.6-ext using patch by David Bogar - Penn State, ARL

    Added capability to analyse both pumps and turbines using a flag in the
    controlDict by Bryan Lewis - Penn State

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    turboPerformance.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_turboPerformance_H
#define functionObjects_turboPerformance_H

#include "turboPerformance/fluidPower.H"
#include "rhoThermo/rhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class turboPerformance Declaration
\*---------------------------------------------------------------------------*/

class turboPerformance
:
    public fluidPower
{
public:

    //- Machine types
        enum machineType
        {
            pump,
            hydraulicReactionTurbine,
            fan,
            compressor,
            gasTurbine
        };

        //- Machine types
        static const NamedEnum<machineType, 5> machineTypeNames_;

    //- Machine types
        enum weightingType
        {
            area,
            flux
        };

        //- Machine types
        static const NamedEnum<weightingType, 2> weightingTypeNames_;

private:

    // Private data

        //- Machine types
        machineType machineType_;

        //- Weighting types
        weightingType weightingType_;

        //- moving wall sections
        const label movingWallSectors_;

        //- Switch to use generic formula for efficiency at compressor and
        // gasTurbine machines. This is valid for every compressible material
        const Switch realGas_;

        //- Switch to print intermidiate quantities for debugging
        //  Useful when the results don't make sense
        const Switch debugInfo_;

        //- Add sector force and moment to 360
        void transformSectorQuantitiesTo360(vector&, vector&) const;

        //- Machine specific output header information
        void writeFileHeaderPumpOrTurbine(Ostream& os) const;
        void writeFileHeaderFan(Ostream& os) const;
        void writeFileHeaderCompGasTurbine(Ostream& os) const;

        //- S from thermo.S(p,T)
        scalar entropy
        (
            fluidThermo&,
            const scalar&,
            const scalar&
        );

        //- h from thermo.HE(p,T)
        scalar enthalpy
        (
            fluidThermo&,
            const scalar&,
            const scalar&
        );

        //- 2nd order FD to find ds/dT for N-R
        scalar dSdT
        (
            fluidThermo&,
            const scalar&,
            const scalar&,
            const scalar&
        );

        //- Calulate stagnation T0 p0 by solving 2x2 system
        //  s1 = s2  and h1 = h2
        //  Using FD to find ds/dt ds/dp dh/dt dh/dp
        //  N-R 2x2 to converge and find T0 p0
        void calculateStagnation
        (
            fluidThermo& thermo,
            const scalar& p,
            const scalar& t,
            const scalar& ht,
            const scalar& s,
            scalar& p0,
            scalar& T0
        );

        //- Calulate Tisentropic based on s1 = s2(Tis, p2) = const
        //  Using Newton-Raphson and FD to find derivative
        void calculateTis
        (
            fluidThermo& thermo,
            const scalar& p2,
            const scalar& s1,
            const label& maxIter,
            const scalar& tolerance,
            scalar& T2ist
        );

        //- Machine specific
        void executePumpOrTurbine
        (
            const vector&,
            const vector&,
            const vector&,
            const scalar&
        );
        void executeFan
        (
            const vector&,
            const vector&,
            const vector&,
            const scalar&
        );

        //- Using ideal gas relations
        void executeCompGasTurbine
        (
            const vector&,
            const vector&,
            const vector&,
            const scalar&
        );

        //- Generic - based on the material library
        void executeCompGasTurbineRealGas
        (
            const vector&,
            const vector&,
            const vector&,
            const scalar&
        );


protected:

    //- Output file header information
    virtual void writeFileHeader(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("turboPerformance");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        turboPerformance
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        turboPerformance(const turboPerformance&) = delete;


    //- Destructor
    virtual ~turboPerformance();


    // Member Functions

        //- Execute, currently does nothing
        virtual bool execute();

        //- Read the forces data
        virtual bool read(const dictionary&);

        //- Calculate the forces and write
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const turboPerformance&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
