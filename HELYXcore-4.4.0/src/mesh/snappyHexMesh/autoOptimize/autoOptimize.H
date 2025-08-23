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
    (c) 2016 Engys Ltd.

Class
    Foam::autoOptimize

Description
    Perform optimisation of mesh

SourceFiles
    autoOptimize.C

\*---------------------------------------------------------------------------*/

#ifndef autoOptimize_H
#define autoOptimize_H

#include "fvMesh/fvMesh.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class autoOptimize Declaration
\*---------------------------------------------------------------------------*/

class autoOptimize
{

public:

        //- Runtime type information
        TypeName("autoOptimize");

        // Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            autoOptimize,
            dictionary,
            (
                fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );

        // Constructors

        autoOptimize
        (
            fvMesh& mesh,
            const dictionary& dict
        );

        //- Destructor
        virtual ~autoOptimize();

        // Selectors

        //- Select null constructed
        static autoPtr<autoOptimize> New
        (
            fvMesh& mesh,
            const dictionary& dict
        );

        //- Update input points field
        virtual void movePoints(pointField& newPoints) = 0;

        //- optimise mesh points
        virtual void optimize() = 0;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
