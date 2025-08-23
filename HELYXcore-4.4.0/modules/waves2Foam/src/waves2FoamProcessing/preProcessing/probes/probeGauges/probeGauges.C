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

#include "probeGauges.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(probeGauges, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void probeGauges::writeVTKFormat
(
    const word& name,
    const pointField& pp
)
{
    autoPtr<OFstream> vtk;

    // Writing the lines as VTK format
    vtk.reset(new OFstream("waveGaugesNProbes/" + name + ".vtk"));

    // Writing header
    vtk() << "# vtk DataFile Version 3.0" << nl << "vtk output" << nl
          << "ASCII" << nl << "DATASET POLYDATA" << endl;

    // Writing points
    vtk() << "POINTS " << pp.size() << " float" << endl;

    forAll(pp, pointi)
    {
        point p( pp[pointi] );
        vtk() << p.x() << " " << p.y() << " " << p.z() << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


probeGauges::probeGauges
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),

    gaugeDict_(dict)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void probeGauges::evaluate(const word& name)
{
	// Evaluate the point distribution
	autoPtr<Foam::pointDistributions> pd
	(
		Foam::pointDistributions::New(mesh_, gaugeDict_)
	);

	pointField pp(pd->evaluate());

	// Get the additional information for the sampling
	label interval(gaugeDict_.lookup<label>("outputInterval"));

	wordList fields(gaugeDict_.lookup("fields"));

	// Write the file to be included in controlDict
    autoPtr<OFstream> gauges;

    gauges.reset(new OFstream( "waveGaugesNProbes/" + name + "_controlDict" ));

    gauges() << incrIndent << indent << name << nl << indent
             << token::BEGIN_BLOCK << nl << incrIndent;
    gauges() << indent << "type               probes;" << nl;
    gauges() << indent << "functionObjectLibs ( \"libsampling.so\" );" << nl;
    gauges() << indent << "enabled            true;" << nl;
    gauges() << nl;
    gauges() << indent << "outputControl      timeStep;" << nl;
    gauges() << indent << "outputInterval     " << interval << ";" << nl << nl;
    gauges() << nl;
    gauges() << indent << "probeLocations" << nl << indent << "(" << nl;
    forAll(pp, pointi)
    {
        gauges() << indent << "    " << pp[pointi] << nl;
    }
    gauges() << indent << ");" << nl;
    gauges() << nl;
    gauges() << indent << "fields" << nl << indent << "(" << nl;
    forAll(fields, fieldi)
    {
        gauges() << indent << "    " << fields[fieldi] << nl;
    }
    gauges() << indent << ");" << endl;

    gauges() << decrIndent << indent << token::END_BLOCK << decrIndent << nl;
    gauges() << nl;

    // Write the point locations in VTK-format
    if (gaugeDict_.lookupOrDefault<Switch>("writeVTK", true))
    {
        writeVTKFormat(name, pp);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
