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
    (c) 2015-2021 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::regionModels::surfaceFilmModels::patchEjection

Description
    Remove and eject the mass in the film as it passes over the selected
    patches.

SourceFiles
    patchEjection.C

\*---------------------------------------------------------------------------*/

#ifndef patchEjection_H
#define patchEjection_H

#include "surfaceFilmModels/submodels/kinematic/ejectionModel/ejectionModel/ejectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class patchEjection Declaration
\*---------------------------------------------------------------------------*/

class patchEjection
:
    public ejectionModel
{
protected:

        //- Stable film thickness - mass only removed if thickness execeeds
        //  this threhold value
        scalar deltaStable_;

        //- List of patch IDs at which the film is removed
        labelList patchIDs_;

        //- Ejected mass for each patch at which the film is removed
        scalarField patchEjectedMasses_;


public:

    //- Runtime type information
    TypeName("patchEjection");


    // Constructors

        //- Construct from surface film model
        patchEjection(surfaceFilmRegionModel& film, const dictionary& dict);

        //- Disallow default bitwise copy construction
        patchEjection(const patchEjection&) = delete;


    //- Destructor
    virtual ~patchEjection();


    // Member Functions

        //- Correct
        virtual void correct
        (
            scalarField& availableMass,
            scalarField& massToEject,
            scalarField& diameterToEject
        );

        //- Accumulate the total mass ejected for the patches into the
        //  scalarField provided
        virtual void patchEjectedMassTotals
        (
            scalarField& patchMasses
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const patchEjection&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
