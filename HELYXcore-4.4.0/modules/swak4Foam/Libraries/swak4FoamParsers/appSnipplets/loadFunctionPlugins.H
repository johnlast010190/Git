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
    (c) 1991-2005 OpenCFD Ltd. ICE Stroemungsfoschungs GmbH

Description
    Load function plugins if specified on the command line

Contributors/Copyright:
    2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

if (args.options().found("functionPlugins")) {
    Info<< "\nLoading additional function plugins" << endl;

    string theOption(args.options()["functionPlugins"]);
    std::string::size_type start=0;
    std::string::size_type end=0;

    while (start<theOption.length()) {
        end=theOption.find(',',start);
        if (end==std::string::npos) {
            end=theOption.length();
        }
        word name(
            theOption.substr(
                start,end-start
            )
        );
        start=end+1;

        fileName theLib("libswak"+name+"FunctionPlugin.so");
        Info<< "  Loading " << theLib << endl;

#ifdef FOAM_DLLIBRARY_USES_STATIC_METHODS
        dlLibraryTable::open(theLib);
#else
        table.open(theLib);
#endif
    }

}

// ************************************************************************* //
