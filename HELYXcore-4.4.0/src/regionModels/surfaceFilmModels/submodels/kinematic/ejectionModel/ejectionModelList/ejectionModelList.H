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
    Foam::regionModels::surfaceFilmModels::ejectionModelList

Description
    List container for film ejection models

SourceFiles
    ejectionModelList.C

\*---------------------------------------------------------------------------*/

#ifndef ejectionModelList_H
#define ejectionModelList_H

#include "containers/Lists/PtrList/PtrList.H"
#include "surfaceFilmModels/submodels/kinematic/ejectionModel/ejectionModel/ejectionModel.H"
#include "surfaceFilmModels/submodels/filmSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                    Class ejectionModelList Declaration
\*---------------------------------------------------------------------------*/

class ejectionModelList
:
    public PtrList<ejectionModel>,
    public filmSubModelBase
{
private:

    // Private data

        //- List of mass ejected per patch
        scalarField massEjected_;


public:

    // Constructors

        //- Construct null
        ejectionModelList(surfaceFilmRegionModel& film);

        //- Construct from type name, dictionary and surface film model
        ejectionModelList
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        ejectionModelList(const ejectionModelList&) = delete;


    //- Destructor
    virtual ~ejectionModelList();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct
            (
                scalarField& availableMass,
                volScalarField& massToEject,
                volScalarField& diameterToEject
            );


        // I-O

            //- Provide some info
            virtual void info(Ostream& os);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ejectionModelList&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
