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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::inversePointDistanceDiffusivity

Description
    Inverse distance to the given patches motion diffusivity.

SourceFiles
    inversePointDistanceDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef inversePointDistanceDiffusivity_H
#define inversePointDistanceDiffusivity_H

#include "motionDiffusivity/uniform/uniformDiffusivity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class inversePointDistanceDiffusivity Declaration
\*---------------------------------------------------------------------------*/

class inversePointDistanceDiffusivity
:
    public uniformDiffusivity
{
    // Private data

        //- Patches selected to base the distance on
        wordReList patchNames_;


public:

    //- Runtime type information
    TypeName("inversePointDistance");


    // Constructors

        //- Construct for the given fvMesh and data Istream
        inversePointDistanceDiffusivity(const fvMesh& mesh, Istream& mdData);

        //- Disallow default bitwise copy construction
        inversePointDistanceDiffusivity
        (
            const inversePointDistanceDiffusivity&
        ) = delete;


    //- Destructor
    virtual ~inversePointDistanceDiffusivity();


    // Member Functions

        //- Correct the motion diffusivity
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const inversePointDistanceDiffusivity&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
