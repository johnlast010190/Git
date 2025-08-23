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
    Foam::transferModelList

Description
    List container for film transfer models

SourceFiles
    transferModelList.C

\*---------------------------------------------------------------------------*/

#ifndef transferModelList_H
#define transferModelList_H

#include "containers/Lists/PtrList/PtrList.H"
#include "surfaceFilmModels/submodels/kinematic/transferModels/transferModel/transferModel.H"
#include "surfaceFilmModels/submodels/filmSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                    Class transferModelList Declaration
\*---------------------------------------------------------------------------*/

class transferModelList
:
    public PtrList<transferModel>,
    public filmSubModelBase
{
    // Private data

        //- List of mass transferred per patch
        scalarField massTransferred_;


public:

    // Constructors

        //- Construct null
        transferModelList(surfaceFilmRegionModel& film);

        //- Construct from type name, dictionary and surface film model
        transferModelList
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construct
        transferModelList(const transferModelList&) = delete;


    //- Destructor
    virtual ~transferModelList();


    // Member Functions

        //- Correct kinematic transfers
        virtual void correct
        (
            scalarField& availableMass,
            volScalarField& massToTransfer,
            volVectorField& momentumToTransfer
        );

        //- Correct kinematic and thermodynamic transfers
        virtual void correct
        (
            scalarField& availableMass,
            volScalarField& massToTransfer,
            volVectorField& momentumToTransfer,
            volScalarField& energyToTransfer
        );

        //- Provide some info
        virtual void info(Ostream& os);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const transferModelList&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
