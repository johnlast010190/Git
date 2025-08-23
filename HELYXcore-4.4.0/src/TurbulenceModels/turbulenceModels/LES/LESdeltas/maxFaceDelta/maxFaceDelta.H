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
    (c) 2011 OpenFOAM Foundation
    (c) 2010-2024 Engys Ltd.

Class
    Foam::maxFaceDelta

Description
    maxFaceDelta returns the maximum mean thickness of the cell by comparing
    the mean thickness in the direction of each face.
    max_i(2*V/sum_j(mag((Cf - C) & Sf_j))).
    Valid for arbitrary cell shapes.


SourceFiles
    maxFaceDelta.C

\*---------------------------------------------------------------------------*/

#ifndef maxFaceDeltaDelta_H
#define maxFaceDeltaDelta_H

#include "LES/LESdeltas/LESdelta/LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class maxFaceDelta Declaration
\*---------------------------------------------------------------------------*/

class maxFaceDelta
:
    public LESdelta
{
    // Private data

        scalar deltaCoeff_;

        //- recompute deltas
        bool recalcDelta_;


    // Private Member Functions

        // Calculate the delta values
        void calcDelta();


public:

    //- Runtime type information
    TypeName("maxFaceDelta");


    // Constructors

        //- Construct from name, mesh and IOdictionary
        maxFaceDelta
        (
            const word& name,
            const turbulenceModel& turbulence,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        maxFaceDelta(const maxFaceDelta&) = delete;


    //- Destructor
    virtual ~maxFaceDelta()
    {}


    // Member Functions

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&);

        // Correct values
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const maxFaceDelta&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
