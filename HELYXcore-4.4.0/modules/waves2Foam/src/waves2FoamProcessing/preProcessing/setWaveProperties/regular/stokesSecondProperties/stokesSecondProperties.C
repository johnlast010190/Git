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
    (c) held by original author

\*---------------------------------------------------------------------------*/

#include "stokesSecondProperties.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(stokesSecondProperties, 0);
addToRunTimeSelectionTable
(
    setWaveProperties,
    stokesSecondProperties,
    setWaveProperties
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


stokesSecondProperties::stokesSecondProperties
(
    const Time& rT,
    dictionary& dict,
    bool write
)
:
    setWaveProperties(rT, dict, write),
    sfp_( rT, dict, false, "")
{
    Info<< "\nConstructing: " << this->type() << endl;

    period_ = dict.lookup<scalar>("period");
    depth_  = dict.lookup<scalar>("depth");
    omega_  = 2.0*PI_/period_ ;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void stokesSecondProperties::set(Ostream& os)
{
    scalar k = sfp_.linearWaveNumber();

    // Write the beginning of the sub-dictionary
    writeBeginning( os );

    // Write the already given parameters
    writeGiven( os, "waveType" );

    if (dict_.found("Tsoft"))
    {
        writeGiven( os, "Tsoft");
    }

    writeGiven( os, "depth" );
    writeGiven( os, "period" );
    writeGiven( os, "direction" );
    writeGiven( os, "phi");
    writeGiven( os, "height");

    if (write_)
    {
        vector direction(dict_.lookup<vector>("direction"));
        direction /= Foam::mag(direction);
        direction *= k;

        writeDerived(os, "waveNumber", direction);
        writeDerived(os, "omega", sfp_.omega());
    }

    writeGiven( os, "debug");

    // Write the relaxation zone
    writeRelaxationZone( os );

    // Write the closing bracket
    writeEnding( os );

    scalar H = dict_.lookup<scalar>("height");
    scalar h = dict_.lookup<scalar>("depth");

    scalar a1 = H/2.0;
    scalar a2 = 1.0/16.0*k * sqr(H)
        *(3.0/Foam::pow(Foam::tanh(k*h),3.0) - 1.0/Foam::tanh(k*h));

    if (Switch(dict_.lookup("debug")))
    {
        Info<< nl << "The wave amplitudes are:\n" << tab << "  a1 = "
             << tab << a1
             << nl << tab << "  a2 = " << tab << a2
             << nl << tab << "4 a2 = " << tab << 4.0*a2
             << " (Validity criterion) " << endl;
    }

    if (a1 < 4.0*a2)
    {
        Info<< a1 << tab << 4.0*a2 << endl;

        WarningIn
        (
            "void stokesSecondProperties::set(Ostream& os)"
        ) << endl << "The validity of Stokes second order is violated." << endl
          << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
