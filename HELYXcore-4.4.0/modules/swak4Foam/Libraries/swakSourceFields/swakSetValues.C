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
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "include/swak.H"

#ifdef FOAM_HAS_BASICSOURCES

#include "SwakSetValue.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

#if defined(FOAM_HAS_FVOPTIONS) && !defined(FOAM_FVOPTION_MAKE_NOT_IN_NAMESPACE)
    namespace fv {
#endif

        makeSwakFvOption(SwakSetValue, scalar);
        makeSwakFvOption(SwakSetValue, vector);
        makeSwakFvOption(SwakSetValue, sphericalTensor);
        makeSwakFvOption(SwakSetValue, symmTensor);
        makeSwakFvOption(SwakSetValue, tensor);

#if defined(FOAM_HAS_FVOPTIONS) && !defined(FOAM_FVOPTION_MAKE_NOT_IN_NAMESPACE)
    }
#endif

}

#endif

// ************************************************************************* //
