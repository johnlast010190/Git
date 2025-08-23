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


Description
    Define types that were renamed in OpenFOAM 2.2

SourceFiles

Contributors/Copyright:
    2012-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef SwakFvOptionMacroHeader_H
#define SwakFvOptionMacroHeader_H

#include "swak.H"

#ifdef FOAM_HAS_FVOPTIONS

#include "cfdTools/general/fvOptions/fvOption.H"
#include "cfdTools/general/fvOptions/makeFvOption.H"

#if (FOAM_VERSION4SWAK >= VERSION_NR(3,0,0)) && !defined(FOAM_VERSION4SWAK_IS_EXTEND)

#include "cellSetOption/cellSetOption.H"

namespace Foam {
    typedef fv::cellSetOption swakFvOptionType;
}

#define FOAM_FVOPTION_MAKE_NOT_IN_NAMESPACE

#else



namespace Foam {
    typedef fv::option swakFvOptionType;
}

#endif

#define makeSwakFvOption makeFvOption

#if (FOAM_VERSION4SWAK >= VERSION_NR(2,3,1)) && !defined(FOAM_VERSION4SWAK_IS_EXTEND)
#define FOAM_FVOPTION_HAS_ADDITIONAL_ADDSUP
#endif

#else

#include "basicSource.H"
#include "makeBasicSource.H"

namespace Foam {
    typedef basicSource swakFvOptionType;
}

#define makeSwakFvOption makeBasicSource

#endif

#endif

// ************************************************************************* //
