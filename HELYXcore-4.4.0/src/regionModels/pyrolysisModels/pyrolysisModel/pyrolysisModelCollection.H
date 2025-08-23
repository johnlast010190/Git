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
    (c) 2011-2014 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::regionModels::pyrolysisModelCollection

Description
    A centralized pyrolysis collection.

    Container class for a set of pyrolysis with functions implemented
    to loop over the functions for each type.

SourceFiles
    pyrolysisModelCollection.C

\*---------------------------------------------------------------------------*/

#ifndef pyrolysisModelCollection_H
#define pyrolysisModelCollection_H

#include "containers/Lists/PtrList/PtrList.H"
#include "pyrolysisModels/pyrolysisModel/pyrolysisModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class fvMesh;

namespace regionModels
{

/*---------------------------------------------------------------------------*\
                  Class pyrolysisModelCollection Declaration
\*---------------------------------------------------------------------------*/

class pyrolysisModelCollection
:
    public PtrList<pyrolysisModel>
{
public:

    // Runtime type information
    TypeName("pyrolysisModelCollection");


    // Constructors

        //- Construct from mesh
        pyrolysisModelCollection(const fvMesh&);

        //- Disallow default bitwise copy construction
        pyrolysisModelCollection(const pyrolysisModelCollection&) = delete;


    //- Destructor
    virtual ~pyrolysisModelCollection();


    // Member Functions

        //- Evolve the pyrolysis equation regions
        virtual void evolveRegion();

        //- Evolve regions
        virtual void evolve();

        //- Provide some feedback from pyrolysis regions
        virtual void info();

        //- Return max diffusivity allowed in the solid
        virtual scalar maxDiff() const;

        //- Mean diffusion number of the solid regions
        virtual scalar solidRegionDiffNo() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pyrolysisModelCollection&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam


#endif

// ************************************************************************* //
