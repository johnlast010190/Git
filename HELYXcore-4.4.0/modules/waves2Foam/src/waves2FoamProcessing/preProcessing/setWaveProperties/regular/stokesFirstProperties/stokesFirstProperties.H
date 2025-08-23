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
    Foam::stokesFirstProperties

Description

    A description of the general wave framework is given in

    @article { jacobsenFuhrmanFredsoe2011,
        Author = {Jacobsen, N G and Fuhrman, D R and Freds\o{}e, J},
        title = {{A Wave Generation Toolbox for the Open-Source CFD Library: OpenFoam\textregistered{}}},
        Journal = {{Int. J. for Numer. Meth. Fluids}},
        Year = {2012},
        Volume = {70},
        Number = {9},
        Pages = {1073-1088},
        DOI = {{10.1002/fld.2726}},
    }

SourceFiles
    setWaveProperties.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef stokesFirstProperties_H
#define stokesFirstProperties_H


#include "preProcessing/setWaveProperties/setWaveProperties.H"

namespace Foam
{

class stokesFirstProperties
:
    public setWaveProperties
{
private:

    scalar period_;
    scalar omega_;
    scalar depth_;

public:

        //- Runtime type information
        TypeName("stokesFirstProperties");

        stokesFirstProperties
        (
            const Time&,
            dictionary&
        );

        stokesFirstProperties
        (
            const Time&,
            dictionary&,
            bool
        );

        stokesFirstProperties
        (
            const Time&,
            dictionary&,
            bool,
            word
        );

        // Method

        virtual void set( Ostream& );

        scalar linearWaveNumber() const;

        scalar linearWaveNumber(const scalar&, const scalar&);

        const scalar& omega() const { return omega_; };
};

}

#endif

// ************************************************************************* //
