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
    Foam::JONSWAP

Description

SourceFiles
    JONSWAP.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef JONSWAP_H
#define JONSWAP_H

#include "preProcessing/setWaveProperties/regular/stokesFirstProperties/stokesFirstProperties.H"
#include "preProcessing/setWaveProperties/irregular/waveSpectra/waveSpectra.H"

namespace Foam
{

class JONSWAP
:
    public waveSpectra
{
private:

        scalarField spectralValue
        (
            const scalar&,
            const scalar&,
            const scalar&,
            const scalarField&
        ) const;

public:

        //- Runtime type information
        TypeName("JONSWAP");

        JONSWAP
        (
            const Time&,
            dictionary&,
            scalarField&,
            scalarField&,
            scalarField&,
            vectorField&
        );

        virtual void set(Ostream&);

        virtual wordList list();
};

}

#endif

// ************************************************************************* //
