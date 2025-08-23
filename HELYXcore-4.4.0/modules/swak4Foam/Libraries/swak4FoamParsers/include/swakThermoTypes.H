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

#ifndef SwakThermoMacroHeader_H
#define SwakThermoMacroHeader_H

#include "swak.H"

#ifdef FOAM_HAS_FLUIDTHERMO

#include "rhoThermo/rhoThermo.H"

namespace Foam
{
    typedef fluidThermo swakFluidThermoType;
    typedef rhoThermo swakRhoThermoType;
}

#else

#include "basicRhoThermo.H"

namespace Foam {
    typedef basicThermo swakFluidThermoType;
    typedef basicRhoThermo swakRhoThermoType;
}

#endif

#ifndef FOAM_VERSION4SWAK_IS_EXTEND
#define FOAM_RADIATION_HAS_SOOTMODEL
#define FOAM_RADIATION_HAS_PUBLIC_ABSORPTIONMODEL
#else
#define FOAM_BASIC_THERMO_HAS_NO_NEW
#endif

#if (OPENFOAM_PLUS >= 1706) \
    || (defined(FOAM_VERSION4SWAK_IS_OPENFOAM) && FOAM_VERSION4SWAK>=VERSION_NR(5,0,0))
#define FOAM_CHEMISTRYMODEL_HAS_NO_SOURCE_TERM
#endif

#include "helpers/DebugOStream.H"

#endif

// ************************************************************************* //
