/*---------------------------------------------------------------------------* \
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
    (c) 2018-2025 Engys Ltd.
    (c) D. Segersson
    (c) 2011-2017 OpenFOAM Foundation

Class
    Foam::fv::porousAbsorptionEmissionSource

Group
    grpFvOptionsSources

Description

    Local absorption emission model for participating fvDOM radiation model,
    which supports specification on cell porous zones.

    Absorbed fraction along path ds: a_*ILabda(omega_)
    Emitted fraction along path ds: a_*sigma/pi*T^4

    Interaction with the temperature and energy equation is done via
    source terms (see Ru and Rp).

    See also radiativeIntensityRay.C
             fvDOM.C
             radiationModel.C
    for details about equations.

    source
    {
        type            porousAbsorptionEmissionSource;
        active          yes;
        porousAbsorptionEmissionSourceCoeffs
        {
            selectionMode   cellZone;
            cellZone        zone;
            fields          (ILambda_0_0 ILambda_1_0 ILambda_2_0 ILambda_3_0 ILambda_4_0
            ILambda_5_0 ILambda_6_0 ILambda_7_0 T );
            emissivity      0.8;
            E               0.0;
            rhoSolid        2000;
            CpSolid         1000;
            hCoefficient    23;
            solarHeat       false; // if true, the cellZone heats up because of the solar rad
        }
    }

SourceFiles
    porousAbsorptionEmissionSource.C

\*---------------------------------------------------------------------------*/
#ifndef porousAbsorptionEmissionSource_H
#define porousAbsorptionEmissionSource_H

#include "cellSetOption/cellSetOption.H"
#include "sources/derived/absorptionEmissionSource/absorptionEmissionSource.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
              Class porousAbsorptionEmissionSource Declaration
\*---------------------------------------------------------------------------*/

class porousAbsorptionEmissionSource
:
    public absorptionEmissionSource
{
    // Private Data

        //- Temperature of the solid part of the porous medium
        volScalarField Tsolid_;

        //- density of the solid part of the porous medium
        scalar rhoSolid_;

        //- Heat capacity of the solid part of the porous medium
        scalar CpSolid_;

        //- heat transfer coefficient times area: h*A [W/(Km^3)]
        //- alternatively express A in terms of porosity Vfluid/Vtotal
        scalar hCoeff_;


    // Private Member Functions

        //- update Tsolid of porous medium
        void calculateTsolid();

public:

    //- Runtime type information
    TypeName("porousAbsorptionEmissionSource");


    // Constructors

        // - Construct from components
        porousAbsorptionEmissionSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~porousAbsorptionEmissionSource()
    {}

    // Member Functions

        //- Enthalpy equation sources
        void addSupHe(fvMatrix<scalar>& eqn) override;

        //- ILambda equations absorption and emission terms
        void addSupILambda(fvMatrix<scalar>& eqn) override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const porousAbsorptionEmissionSource&) = delete;
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
