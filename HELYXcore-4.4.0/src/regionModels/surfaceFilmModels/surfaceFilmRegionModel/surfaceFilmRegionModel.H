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

Class
    Foam::regionModels::surfaceFilmModels::surfaceFilmRegionModel

Description
    Base class for surface film models

SourceFiles
    surfaceFilmRegionModel.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceFilmRegionModel_H
#define surfaceFilmRegionModel_H

#include "surfaceFilmModels/surfaceFilmModel/surfaceFilmModel.H"
#include "regionModel/singleLayerRegion/singleLayerRegion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class surfaceFilmRegionModel Declaration
\*---------------------------------------------------------------------------*/

class surfaceFilmRegionModel
:
    public surfaceFilmModel,
    public singleLayerRegion
{
    // Private data

        //- Acceleration due to gravity [m/s2]
        const dimensionedVector& g_;


protected:

    // Protected member functions

        //- Read control parameters from dictionary
        virtual bool read();


public:

    //- Runtime type information
    TypeName("surfaceFilmRegionModel");


    // Constructors

        //- Construct from type name, mesh and gravity vector
        surfaceFilmRegionModel
        (
            const word& modelType,
            const fvMesh& mesh,
            const dimensionedVector& g,
            const word& regionType
        );

        //- Disallow default bitwise copy construct
        surfaceFilmRegionModel(const surfaceFilmRegionModel&) = delete;


    //- Destructor
    virtual ~surfaceFilmRegionModel();


    // Member Functions

        // Access

            //- Return the accleration due to gravity
            inline const dimensionedVector& g() const
            {
                return g_;
            }

            //- External hook to add sources to the film
            virtual void addSources
            (
                const label patchi,
                const label facei,
                const scalar massSource,
                const vector& momentumSource,
                const scalar pressureSource,
                const scalar energySource
            ) = 0;


        // Fields

            //- Return the film thickness [m]
            virtual const volScalarField& delta() const = 0;

            //- Return the film coverage, 1 = covered, 0 = uncovered []
            virtual const volScalarField& coverage() const = 0;

            //- Return the film velocity [m/s]
            virtual const volVectorField& U() const = 0;

            //- Return the film surface velocity [m/s]
            virtual tmp<volVectorField::Internal> Us() const = 0;

            //- Return the film wall velocity [m/s]
            virtual tmp<volVectorField::Internal> Uw() const = 0;

            //- Return the film density [kg/m3]
            virtual const volScalarField& rho() const = 0;

            //- Return the film surface tension [N/m]
            virtual tmp<volScalarField> sigma() const = 0;


        // Transfer fields - to the primary region

            //- Return mass transfer source - Eulerian phase only
            virtual tmp<volScalarField> primaryMassTrans() const = 0;

            //- Return the film mass available for transfer
            virtual const volScalarField& cloudMassTrans() const = 0;

            //- Return the parcel diameters originating from film
            virtual const volScalarField& cloudDiameterTrans() const = 0;


        // Evolution

            //- Main driver routing to evolve the region - calls other evolves
            virtual void evolve();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const surfaceFilmRegionModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //