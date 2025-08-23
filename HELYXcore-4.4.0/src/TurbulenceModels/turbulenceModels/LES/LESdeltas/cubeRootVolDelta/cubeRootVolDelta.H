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
    (c) 2016 OpenCFD ltd.
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::cubeRootVolDelta

Description
    Simple cube-root of cell volume delta used in LES models.

SourceFiles
    cubeRootVolDelta.C

\*---------------------------------------------------------------------------*/

#ifndef cubeRootVolDelta_H
#define cubeRootVolDelta_H

#include "LES/LESdeltas/LESdelta/LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class cubeRootVolDelta Declaration
\*---------------------------------------------------------------------------*/

class cubeRootVolDelta
:
    public LESdelta
{
    // Private data

        scalar deltaCoeff_;


public:

    //- Runtime type information
    TypeName("cubeRootVol");


    // Constructors

        //- Construct from name, turbulenceModel and dictionary
        cubeRootVolDelta
        (
            const word& name,
            const turbulenceModel& turbulence,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        cubeRootVolDelta(const cubeRootVolDelta&) = delete;


    //- Destructor
    virtual ~cubeRootVolDelta()
    {}


    // Member Functions

        //- Calculate the delta values
        void calcDelta();

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&);

        // Correct values
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cubeRootVolDelta&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
