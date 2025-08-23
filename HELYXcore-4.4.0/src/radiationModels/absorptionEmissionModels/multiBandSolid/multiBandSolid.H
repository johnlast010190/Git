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
    (c) 2015 OpenCFD Ltd.

Class
    Foam::radiationModels::absorptionEmissionModels::multiBandSolid

Description

    multiBandSolid radiation absorption/emission for solids.


SourceFiles
    multiBandSolid.C

\*---------------------------------------------------------------------------*/

#ifndef multiBandSolid_H
#define multiBandSolid_H

#include "absorptionEmissionModels/absorptionEmissionModel/absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{

/*---------------------------------------------------------------------------*\
                  Class multiBandSolid Declaration
\*---------------------------------------------------------------------------*/

class multiBandSolid
:
    public absorptionEmissionModel
{
public:

    // Public data

        //- Maximum number of bands
        static const int maxBands_ = 5;


private:

    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //-  Absorption coefficients
        scalarList absCoeffs_;

        //-  Emissivity coefficients
        scalarList emiCoeffs_;

        //- Bands
        label nBands_;


public:

    //- Runtime type information
    TypeName("multiBandSolid");


    // Constructors

        //- Construct from components
        multiBandSolid
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& typeNameDerived = typeName
        );


    //- Destructor
    virtual ~multiBandSolid();


    // Member Functions

        //- Absorption coefficient
        tmp<volScalarField> aCont(const label bandI) const;

        //- Emission coefficient
        tmp<volScalarField> eCont(const label bandI) const;

        //- Emission contribution
        tmp<volScalarField> ECont(const label bandI) const;


        inline bool isGrey() const
        {
            return false;
        }

        //- Number of bands
        inline label nBands() const
        {
            return nBands_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace absorptionEmissionModels
} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
