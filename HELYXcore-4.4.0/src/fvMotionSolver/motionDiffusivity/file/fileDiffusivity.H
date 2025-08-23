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
    Foam::fileDiffusivity

Description
    Motion diffusivity read from given file name.

SourceFiles
    fileDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef fileDiffusivity_H
#define fileDiffusivity_H

#include "motionDiffusivity/motionDiffusivity/motionDiffusivity.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fileDiffusivity Declaration
\*---------------------------------------------------------------------------*/

class fileDiffusivity
:
    public motionDiffusivity
{

protected:

    // Protected data

        surfaceScalarField faceDiffusivity_;


public:

    //- Runtime type information
    TypeName("file");


    // Constructors

        //- Construct for the given fvMesh and data Istream
        fileDiffusivity(const fvMesh& mesh, Istream& mdData);

        //- Disallow default bitwise copy construction
        fileDiffusivity(const fileDiffusivity&) = delete;


    //- Destructor
    virtual ~fileDiffusivity();


    // Member Functions

        //- Return diffusivity field
        virtual tmp<surfaceScalarField> operator()() const
        {
            return faceDiffusivity_;
        }

        //- Do not correct the motion diffusivity
        virtual void correct()
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fileDiffusivity&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
