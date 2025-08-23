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

Application
    probesNGauges

Description
    This utility writes the needed input files for wave gauges to ease the
    definition of these.

    In a later stage it is intended to be extended to be used for write the
    location of point probes as well.

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"

#if EXTBRANCH==1
    #if 310<OFVERSION
        #include "db/Time/foamTime.H"
    #else
        #include "db/Time/Time.H"
    #endif
#elif OFPLUSBRANCH==1
    #include "db/Time/Time.H"
#else
    #include "db/Time/Time.H"
#endif

//#if EXTBRANCH==1 && OFVERSION>310
//    #include "foamTime.H"
//#else
//    #include "Time.H"
//#endif

#include "fvMesh/fvMesh.H"

#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "preProcessing/probes/waveGauges/waveGauges.H"
#include "preProcessing/probes/probeGauges/probeGauges.H"

using namespace Foam;

int main(int argc, char *argv[])
{

#   include "include/addTimeOptions.H"
#   include "include/setRootCase.H"

#   include "include/createTime.H"
#   include "include/createMesh.H"

    // Needed by e.g. the reflection analysis
#   include "cfdTools/general/include/readGravitationalAcceleration.H"
    Info<< "\n";

    mkDir("waveGaugesNProbes");

    IOdictionary probeDefs
    (
        IOobject
        (
            "probeDefinitions",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    wordList toc(probeDefs.toc());

    forAll(toc, item)
    {
        word name(toc[item]);

        if (probeDefs.isDict(name))
        {
            const dictionary& dict(probeDefs.subDict(name));

            if (word(dict.lookup("type")) == "waveGauge")
            {
                waveGauges wg(mesh, dict);

                wg.evaluate(name);
            }
            else if (word(dict.lookup("type")) == "probeGauge")
            {
            	probeGauges pg(mesh, dict);

            	pg.evaluate(name);
            }
            else
            {
                Info<< "Probe-type: '"
                     << word(dict.lookup("type"))
                     << "' not yet implemented" << endl;
            }
        }
    }

    Info<< nl << "End" << endl;

    return 0;
}

