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
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::fv::constantHeatTransfer

Group
    grpFvOptionsSources

Description
    Constant heat transfer model. htcConst [W/m2/K] and area/volume [1/m]
    must be provided.

\*---------------------------------------------------------------------------*/

#ifndef constantHeatTransfer_H
#define constantHeatTransfer_H

#include "sources/interRegion/interRegionHeatTransfer/interRegionHeatTransferModel/interRegionHeatTransferModel.H"
#include "memory/autoPtr/autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class constantHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class constantHeatTransfer
:
    public interRegionHeatTransferModel
{
private:

    // Private data

        //- Constant heat transfer coefficient [W/m2/K]
        autoPtr<volScalarField> htcConst_;

        //- Area per unit volume of heat exchanger [1/m]
        autoPtr<volScalarField> AoV_;


public:

    //- Runtime type information
      TypeName("constantHeatTransfer");


    // Constructors

        //- Construct from dictionary
        constantHeatTransfer
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~constantHeatTransfer();


    // Public Functions

        //- Calculate the heat transfer coefficient
        virtual void calculateHtc(scalar deltaT);


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
