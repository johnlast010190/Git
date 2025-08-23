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
    (c) 2017 OpenCFD Ltd.
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::inverseDistanceDiffusivity

Description
    Inverse distance to the given patches motion diffusivity.

SourceFiles
    inverseDistanceDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef inverseDistanceDiffusivity_H
#define inverseDistanceDiffusivity_H

#include "motionDiffusivity/uniform/uniformDiffusivity.H"
#include "primitives/strings/lists/wordReList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class inverseDistanceDiffusivity Declaration
\*---------------------------------------------------------------------------*/

class inverseDistanceDiffusivity
:
    public uniformDiffusivity
{
    // Private data

        //- Patches selected to base the distance on
        //  These can contain patch names or regular expressions to search for.
        wordReList patchNames_;

        //- Set a buffer region (optional)
        scalar Lbuff_;


public:

    //- Runtime type information
    TypeName("inverseDistance");


    // Constructors

        //- Construct for the given fvMesh and data Istream
        inverseDistanceDiffusivity(const fvMesh& mesh, Istream& mdData);

        //- Disallow default bitwise copy construction
        inverseDistanceDiffusivity(const inverseDistanceDiffusivity&) = delete;


    //- Destructor
    virtual ~inverseDistanceDiffusivity();


    // Member Functions

        //- Correct the motion diffusivity
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const inverseDistanceDiffusivity&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
