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
    (c) 2017-2021 OpenFOAM Foundation

Class
    Foam::regionModels::surfaceFilmModels::transferModel

Description
    Base class for film transfer models, handling mass transfer between the
    film and the continuous phase.

SourceFiles
    transferModel.C
    transferModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef transferModel_H
#define transferModel_H

#include "surfaceFilmModels/submodels/filmSubModelBase.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "fields/Fields/scalarField/scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class transferModel Declaration
\*---------------------------------------------------------------------------*/

class transferModel
:
    public filmSubModelBase
{
    // Private data

        //- Transferred mass
        scalar transferredMass_;


protected:

    // Protected Member Functions

        //- Add to transferred mass
        void addToTransferredMass(const scalar dMass);

        //- Correct
        void correct();


public:

    //- Runtime type information
    TypeName("transferModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            transferModel,
            dictionary,
            (
               surfaceFilmRegionModel& film,
               const dictionary& dict
            ),
            (film, dict)
        );


    // Constructors

        //- Construct for film
        transferModel(surfaceFilmRegionModel& film);

        //- Construct from type name, dictionary and surface film model
        transferModel
        (
            const word& modelType,
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construct
        transferModel(const transferModel&) = delete;


    // Selectors

        //- Return a reference to the selected ejection model
        static autoPtr<transferModel> New
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict,
            const word& modelType
        );


    //- Destructor
    virtual ~transferModel();


    // Member Functions

        //- Correct kinematic transfers
        virtual void correct
        (
            scalarField& availableMass,
            scalarField& massToTransfer,
            vectorField& momentumToTransfer
        ) = 0;

        //- Correct kinematic and thermodynamic transfers
        virtual void correct
        (
            scalarField& availableMass,
            scalarField& massToTransfer,
            vectorField& momentumToTransfer,
            scalarField& energyToTransfer
        ) = 0;

        //- Return the total mass transferred
        virtual scalar transferredMassTotal() const;

        //- Accumulate the total mass transferred for the patches into the
        //  scalarField provided
        virtual void patchTransferredMassTotals(scalarField& patchMasses) const
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const transferModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
