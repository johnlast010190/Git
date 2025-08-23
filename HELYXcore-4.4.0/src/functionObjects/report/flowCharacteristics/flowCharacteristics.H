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
    (c) 2022-2024 Engys Ltd.

Class
    Foam::functionObjects::flowCharacteristics

Group
    grpReportFunctionObjects

Description
    Computes flow characteristics (length, velocity, dT)

    L=min/max (L_vol, L_ext)
    L_vol = pow(volume, 1.0/3.0);
    L_ext = maximum direction of BCs
    U_vol = maximum magnitude of cell values
    U_BC = maximum magnitude of patch faces
    U_dp = sqrt((p_max - p_min)/rhoAverage)

    Based on the above quantities we can compute the following deltaT:
    dt_ = frequency * L/U as aggressive or conservative.

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    flowCharacteristics.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_flowCharacteristics_H
#define functionObjects_flowCharacteristics_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "db/Time/Time.H"
#include "cfdTools/general/include/fvCFD.H"
#include "db/IOstreams/Fstreams/OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class flowCharacteristics Declaration
\*---------------------------------------------------------------------------*/

class flowCharacteristics
:
    public fvMeshFunctionObject,
    public writeFile
{

    // Private data

        //- Field names
        const word UName_;
        const word pName_;
        const word rhoName_;

        //- Frequency
        const scalar f_;


        //- Computed quantities
        scalar totalVol_;

        //- Characteristic length
        scalar lVol_;
        scalar lExt_;

        //- Characteristic velocities
        scalar UVol_;
        scalar UBC_;
        scalar UDP_;

        //- Timesteps
        scalar dtc_;
        scalar dta_;


    // Private Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Calculate quantities
        void calculate();


public:

    //- Runtime type information
    TypeName("flowCharacteristics");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        flowCharacteristics
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        flowCharacteristics(const flowCharacteristics&) = delete;


    //- Destructor
    virtual ~flowCharacteristics();


    // Member Functions

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Read the parameters
        virtual bool read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const flowCharacteristics&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
