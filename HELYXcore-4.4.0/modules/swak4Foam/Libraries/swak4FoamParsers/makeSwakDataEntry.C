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

 ICE Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "include/swak.H"

#ifndef FOAM_HAS_NO_DATAENTRY_CLASS

#include "swakDataEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam {
#ifdef FOAM_DATAENTRY_IS_NOW_FUNCTION1
#define makeDataEntryType makeFunction1Type
#endif

    makeDataEntryType(swakDataEntry,scalar);
    makeDataEntryType(swakDataEntry,vector);
#ifdef FOAM_DATAENTRY_HAS_TENSOR_INSTANCES
    makeDataEntryType(swakDataEntry,tensor);
    makeDataEntryType(swakDataEntry,sphericalTensor);
    makeDataEntryType(swakDataEntry,symmTensor);
#endif
}

#endif

// ************************************************************************* //
