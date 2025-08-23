/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) ICE Stroemungsfoschungs GmbH

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "stdoutFromCommandProvider.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "db/IOstreams/Fstreams/IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(stdoutFromCommandProvider,0);

    // to keep the macro happy
    typedef dynamicFunctionObjectListProxy::dynamicDictionaryProvider dynamicFunctionObjectListProxydynamicDictionaryProvider;

    addToRunTimeSelectionTable
    (
        dynamicFunctionObjectListProxydynamicDictionaryProvider,
        stdoutFromCommandProvider,
        dictionary
    );

    stdoutFromCommandProvider::stdoutFromCommandProvider(
        const dictionary& dict,
        const dynamicFunctionObjectListProxy &owner
    ):
        dynamicFunctionObjectListProxy::dynamicDictionaryProvider(
            dict,
            owner
        ),
        theCommand_(dict.lookup("theCommand"))
    {}


    string stdoutFromCommandProvider::getDictionaryText() {
        string cmd=theCommand_.expand();

        string buffer;
        FILE *output = popen(cmd.c_str(), "r");
        if (!output) {
            FatalErrorIn("stdoutFromCommandProvider::getDictionaryText()")
                << "Problem executing " << cmd
                    << endl
                    << exit(FatalError);

        }
        int c=fgetc(output);
        while (c!=EOF) {
            buffer+=string(char(c));
            c=fgetc(output);
        }
        if (ferror(output)) {
            FatalErrorIn("stdoutFromCommandProvider::getDictionaryText()")
                << "Problem while executing" << cmd
                    << endl
                    << exit(FatalError);

        }
        return buffer;
    }

} // namespace Foam

// ************************************************************************* //
