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
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::NoSurfaceFilm

Group
    grpLagrangianIntermediateSurfaceFilmSubModels

Description
    Place holder for 'none' option

SourceFiles
    NoSurfaceFilm.C

\*---------------------------------------------------------------------------*/

#ifndef NoSurfaceFilm_H
#define NoSurfaceFilm_H

#include "submodels/Kinematic/SurfaceFilmModel/SurfaceFilmModel/SurfaceFilmModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class NoSurfaceFilm Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoSurfaceFilm
:
    public SurfaceFilmModel<CloudType>
{
protected:

    // Protected data

        //- Convenience typedef for parcel type
        typedef typename CloudType::parcelType parcelType;


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from dictionary
        NoSurfaceFilm(const dictionary&, CloudType&);

        //- Construct copy
        NoSurfaceFilm(const NoSurfaceFilm<CloudType>& dm);

        //- Construct and return a clone
        virtual autoPtr<SurfaceFilmModel<CloudType>> clone() const
        {
            return autoPtr<SurfaceFilmModel<CloudType>>
            (
                new NoSurfaceFilm<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoSurfaceFilm();


    // Member Functions

        // Evaluation

            //- Flag to indicate whether model activates the surface film model
            virtual bool active() const;

            //- Transfer parcel from cloud to surface film
            //  Returns true if parcel is to be transferred
            virtual bool transferParcel
            (
                parcelType& p,
                const polyPatch& pp,
                bool& keepParticle
            );

            //- Set parcel properties
            virtual void setParcelProperties
            (
                parcelType& p,
                const label filmCelli
            ) const;


        // I-O

            //- Write surface film info to stream
            virtual void info(Ostream& os);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "submodels/Kinematic/SurfaceFilmModel/NoSurfaceFilm/NoSurfaceFilm.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
