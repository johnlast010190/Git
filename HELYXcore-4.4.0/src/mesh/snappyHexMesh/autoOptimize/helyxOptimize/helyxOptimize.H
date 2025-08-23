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
    (c) 2017 Engys Ltd.

Class
    Foam::helyxOptimize

Description
    Perform optimisation of mesh usign sphericity metric

SourceFiles
    helyxOptimize.C

\*---------------------------------------------------------------------------*/

#ifndef helyxOptimize_H
#define helyxOptimize_H

#include "fvMesh/fvMesh.H"
#include "hessianMeshOptimization/meshOptimization/hessianMeshOptimization.H"
#include "autoOptimize/autoOptimize.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class helyxOptimize Declaration
\*---------------------------------------------------------------------------*/

class helyxOptimize
:
    public autoOptimize
{

private:

      // Private data
        fvMesh& mesh_;

        dictionary coeffsDict_;

        autoPtr<hessianMeshOptimization> optimPtr_;

public:

    //- Runtime type information
    TypeName("helyxOptimize");


        //- Update input points field
        void movePoints(pointField& newPoints);

        //- optimise mesh points
        void optimize();

        // Constructors
        helyxOptimize
        (
            fvMesh& mesh,
            const dictionary& dict
        );

        //- Destructor
        virtual ~helyxOptimize();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
