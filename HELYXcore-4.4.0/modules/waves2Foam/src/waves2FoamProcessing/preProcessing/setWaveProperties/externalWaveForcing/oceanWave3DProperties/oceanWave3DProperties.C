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

#include "oceanWave3DProperties.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(oceanWave3DProperties, 0);
addToRunTimeSelectionTable
(
    setWaveProperties,
    oceanWave3DProperties,
    setWaveProperties
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


oceanWave3DProperties::oceanWave3DProperties
(
    const Time& rT,
    dictionary& dict,
    bool write
)
:
    setWaveProperties(rT, dict, write)
{
    Info<< "\nConstructing: " << this->type() << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void oceanWave3DProperties::set( Ostream& os )
{
    // Write the beginning of the sub-dictionary
    writeBeginning( os );

    // Write the already given parameters
    writeGiven( os, "waveType" );

    os << endl;

    // Write the intervals
    writeGiven(os, "nIntervals");

    scalarField startTimes("startTimes", dict_, dict_.lookup<label>("nIntervals"));
    writeDerived(os, "startTimes", startTimes);

    scalarField endTimes("endTimes", dict_, dict_.lookup<label>("nIntervals"));
    writeDerived(os, "endTimes", endTimes);

    os << endl;

    // Write the ramp information
    writeGiven(os, "rampInterval");

    if (Switch(dict_.lookup("rampInterval")))
    {
    	writeGiven(os, "Tsoft");
    }

    os << endl;

    // Write the information for the mapping
    writeGiven(os, "mappingZone");

    // Set the translation of the OF-mesh
    writeGiven(os, "translateOpenFoamMesh");

    // Write the closing bracket
    writeEnding( os );
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
