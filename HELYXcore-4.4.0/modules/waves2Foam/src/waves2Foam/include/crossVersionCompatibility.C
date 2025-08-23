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

#include "crossVersionCompatibility.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace waves2Foam
{

word pName()
{
#if EXTBRANCH==1
    return "pd";
#elif OFPLUSBRANCH==1
    return "p_rgh";
#else
    #if OFVERSION<170
        return "pd";
    #else
        return "p_rgh";
    #endif
#endif
}


word aName()
{
#if EXTBRANCH==1
    return "alpha1";
#elif OFPLUSBRANCH==1
    return "alpha.water";
#else
    #if OFVERSION<230
        return "alpha1";
    #else
        return "alpha.water";
    #endif
#endif
}


word waterPhase()
{
#if EXTBRANCH==1
    return "phase1";
#elif OFPLUSBRANCH==1
    return "water";
#else
    #if OFVERSION<230
        return "phase1";
    #else
        return "water";
    #endif
#endif
}


word airPhase()
{
#if EXTBRANCH==1
    return "phase2";
#elif OFPLUSBRANCH==1
    return "air";
#else
    #if OFVERSION<230
        return "phase2";
    #else
        return "air";
    #endif
#endif

}


} // End namespace waves2Foam

} // End namespace Foam

// ************************************************************************* //
