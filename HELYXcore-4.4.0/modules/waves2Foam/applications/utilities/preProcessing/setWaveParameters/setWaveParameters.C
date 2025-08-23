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
    setWaveParameters

Description
    This function loops through all of the sub-dictionaries in waveProperties.input
    dictionary in <root case>/constant. The needed wave parameters are computed
    based on a set of required input data.

    Input parameters are given in the file:

        constant/waveProperties.input

    and the derived wave parameters are written to the file

         constant/waveProperties

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.

Additional information
    Implementation published and validated in the following journal article:

    @article { jacobsenFuhrmanFredsoe2011,
        Author = {Jacobsen, N G and Fuhrman, D R and Freds\o{}e, J},
        title = {{A Wave Generation Toolbox for the Open-Source CFD Library: OpenFoam\textregistered{}}},
        Journal = {{Int. J. for Numer. Meth. Fluids}},
        Year = {2012},
        Volume = {70},
        Number = {9},
        Pages = {1073-1088},
        DOI = {{10.1002/fld.2726}},
    }

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "fvMesh/fvMesh.H"
#include "db/IOstreams/IOstreams/IOstream.H"

#include "preProcessing/setWaveProperties/setWaveProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#include "include/setRootCase.H"

#include "include/createTime.H"

    Info<< "\nReading g" << endl;
    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info<< "\nReading waveProperties\n" << endl;

    IOdictionary waveProperties
    (
        IOobject
        (
            "waveProperties.input",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    IOobject wOut
    (
            "waveProperties",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
    );

    // Write waveProperties with the above computed changes
    OFstream os
    (
        wOut.objectPath(),
#if EXTBRANCH==1
        ios_base::out|ios_base::trunc,
#elif OFPLUSBRANCH==1
        // Nothing to be put here
#else
    #if OFVERSION<170
        ios_base::out|ios_base::trunc,
    #endif
#endif
        IOstream::ASCII,
        IOstream::currentVersion,
        IOstream::UNCOMPRESSED
    );

    // Write the OF banner
    wOut.writeBanner( os );

    // Write the file information. Class name is not correct when
    // using wOut.writeHeader( os ); hence manual entries
    os << IOobject::foamFile << nl;
    os << token::BEGIN_BLOCK << incrIndent << nl;
    os << indent << "version" << tab << IOstream::currentVersion
       << token::END_STATEMENT << nl;
    os << indent << "format" << tab << "ascii;" << nl;
    os << indent << "class" << tab << "dictionary;" << nl;
    os << indent << "object" << tab << "waveProperties;" << nl;
    os.endBlock();

    // Write the divider
    wOut.writeDivider( os );
    os << nl;

    /* Loop over all subdicts in waveProperties. For each of them compute the
       wave parameters relevant for that particular wave theory. */
    wordList toc = waveProperties.toc();

    forAll(toc, item)
    {
        // If a sub-dictionary, then compute parameters and write the subdict
        if (waveProperties.isDict(toc[item]))
        {
            dictionary& sd = waveProperties.subDict(toc[item]);

            autoPtr<setWaveProperties> props
                (
                    setWaveProperties::New(runTime, sd, true)
                );

            props->set( os );
        }
        else
        {
            label Nspaces = 20;

            // Read the entry and write to the dummy output file
            ITstream read = waveProperties.lookup(toc[item]);
            os << toc[item] << token::SPACE;

            for (int i=toc[item].size(); i<Nspaces-1; i++)
            {
                os << token::SPACE;
            }

            forAll(read, ri)
            {
                if (ri < read.size() - 1)
                {
                    os << read[ri] << token::SPACE;
                }
                else
                {
                    os << read[ri];
                }
            }

            os << token::END_STATEMENT << nl << endl;

            // Additional level of check, such that the code does not crash at
            // runTime:
            if (toc[item] == "seaLevel")
            {
            	// Read the magnitude of the sea level
            	scalar sL = waveProperties.lookup<scalar>("seaLevel");

            	// If the Switch seaLevelAsReference is not found _and_ the
            	// magnitude of the sea level differs from 0 (zero), stop the
            	// evaluation of the wave parameters
            	if
            	(
            	    !waveProperties.found("seaLevelAsReference") &&
            	    SMALL < Foam::mag(sL)
            	)
            	{
            		// This merely looks up the string, it will not be found
            		// and the user is forced to correct waveProperties.input,
            		// before any execution is possible.
            		waveProperties.lookup("seaLevelAsReference");
            	}
            }
        }
    }

    // Write end divider
    wOut.writeEndDivider(os);

    // End
    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
