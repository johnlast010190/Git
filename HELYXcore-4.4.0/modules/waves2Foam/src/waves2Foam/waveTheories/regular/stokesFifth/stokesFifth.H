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
    Foam::waveTheories::stokesFifth

Description
    Implementation of a fifth order stokes theory for monochromatic surface
    water waves.

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
    stokesFifth.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.

TODO:
    - Compute detadt
    - Compute dpdx

\*---------------------------------------------------------------------------*/

#ifndef stokesFifth_H
#define stokesFifth_H

#include "waveTheories/waveTheory/waveTheory.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

/*---------------------------------------------------------------------------*\
                           Class stokesFifth Declaration
\*---------------------------------------------------------------------------*/

class stokesFifth
:
    public waveTheory
{

protected:

    // Protected data
        scalar H_;
        scalar h_;
        scalar omega_;
        scalar period_;
        scalar phi_;
        vector k_;
        scalar K_;
        scalar Tsoft_;

        // A-coefficients
        scalar A11_;
        scalar A22_;
        scalar A31_;
        scalar A33_;
        scalar A42_;
        scalar A44_;
        scalar A51_;
        scalar A53_;
        scalar A55_;

        // B-coefficients
        scalar B22_;
        scalar B31_;
        scalar B42_;
        scalar B44_;
        scalar B53_;
        scalar B55_;

        // C-coefficients
        scalar C0_;
        scalar C2_;
        scalar C4_;

        // D-coefficients
        scalar D2_;
        scalar D4_;

        // E-coefficients
        scalar E2_;
        scalar E4_;

    // Protected member functions

        scalar factor(const scalar&) const;
public:

    //- Runtime type information
    TypeName("stokesFifth");

    // Constructors

        //- from components
        stokesFifth
        (
            const word&,
            const fvMesh& mesh_
        );


    // Destructor

        ~stokesFifth()
        {}


    // Member Functions
        virtual Switch conflictTSoftInitialise() const
        {
            if (Tsoft_ > 0)
            {
                return true;
            }
            else
            {
                return false;
            }
        };

        void printCoeffs();
        void setCoefficients();

        virtual scalar eta
        (
            const point&,
            const scalar&
        ) const;

        virtual scalar ddxPd
        (
            const point&,
            const scalar&,
            const vector&
        ) const;

        virtual vector U
        (
            const point&,
            const scalar&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
