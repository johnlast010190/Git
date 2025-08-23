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
    (c) 1991-2008 OpenCFD Ltd.

Class
    Foam::SubsetValueExpressionDriver

Description
    print a little banner with the swak-version to ease support

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "printSwakVersion.H"
#include "include/swak.H"

namespace Foam {
    void printSwakVersion() {
        Info<< "swakVersion: " SWAK_VERSION_STRING;
#ifdef SWAK_VERSION_EXTENSION
        Info<< " - " SWAK_VERSION_EXTENSION;
#endif
        Info<< " (Release date: " SWAK_RELEASE_DATE ")";
#ifdef SWAK_HGBRANCH
        Info<< " - HG Branch: " SWAK_HGBRANCH;
#endif
        Info<< endl;
        Info<< "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //";
        Info<< endl;
    }
}
