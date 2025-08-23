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

Class
    Foam::focusingPhase

Description
    Set phases such that the irregular wave train focuses at a given time at
    a given location.

SourceFiles
    focusingPhase.C

Author
    Niels Gjøl Jacobsen, Deltares.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef focusingPhase_H
#define focusingPhase_H

#include "preProcessing/setWaveProperties/irregular/waveSpectra/phases/phases.H"
#include <time.h>

namespace Foam
{

class focusingPhase
:
    public phases
{
private:

        scalar focusTime_;
        vector focusPoint_;

public:

        //- Runtime type information
        TypeName("focusingPhase");

        focusingPhase
        (
            const Time&,
            dictionary&
        );

        virtual scalar phase(const scalar&, const vector&);
};

}

#endif

// ************************************************************************* //
