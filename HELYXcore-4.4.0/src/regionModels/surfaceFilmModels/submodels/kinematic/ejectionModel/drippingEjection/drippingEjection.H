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
    (c) 2011-2021 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::regionModels::surfaceFilmModels::drippingEjection

Description
    Film dripping mass transfer model.

    If the film mass exceeds that needed to generate a valid parcel, the
    equivalent mass is removed from the film.

    New parcel diameters are sampled from a PDF.

SourceFiles
    drippingEjection.C

\*---------------------------------------------------------------------------*/

#ifndef drippingEjection_H
#define drippingEjection_H

#include "surfaceFilmModels/submodels/kinematic/ejectionModel/ejectionModel/ejectionModel.H"
#include "distributionModel/distributionModel.H"
#include "primitives/random/Random/Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                    Class drippingEjection Declaration
\*---------------------------------------------------------------------------*/

class drippingEjection
:
    public ejectionModel
{
protected:

    // Protected data

        //- Stable film thickness - drips only formed if thickness
        //  execeeds this threhold value
        scalar deltaStable_;

        //- Number of particles per parcel
        scalar particlesPerParcel_;

        //- Random number generator
        Random rndGen_;

        //- Parcel size PDF model
        const autoPtr<distributionModel>
            parcelDistribution_;

        //- Diameters of particles to eject into the dripping
        scalarList diameter_;


public:

    //- Runtime type information
    TypeName("drippingEjection");


    // Constructors

        //- Construct from surface film model
        drippingEjection
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        drippingEjection(const drippingEjection&) = delete;


    //- Destructor
    virtual ~drippingEjection();


    // Member Functions

        //- Correct
        virtual void correct
        (
            scalarField& availableMass,
            scalarField& massToEject,
            scalarField& diameterToEject
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const drippingEjection&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
