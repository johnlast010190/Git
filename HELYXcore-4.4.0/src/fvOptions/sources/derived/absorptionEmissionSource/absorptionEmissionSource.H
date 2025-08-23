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
    (c) 2018-2024 Engys Ltd.
    (c) D. Segersson
    (c) 2011-2017 OpenFOAM Foundation

Class
    Foam::fv::absorptionEmissionSource

Group
    grpFvOptionsSources

Description

    Local absorption emission model for participating fvDOM radiation model,
    which supports specification on cell zones.

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
        type            absorptionEmissionSource;
        active          yes;
        absorptionEmissionSourceCoeffs
        {
            selectionMode   cellZone;
            cellZone        zone;
            fields          (ILambda_0_0 ILambda_1_0 ILambda_2_0 ILambda_3_0 ILambda_4_0
            ILambda_5_0 ILambda_6_0 ILambda_7_0 T );
            emissivity      0.8;
            E               0.0;
            solarHeat       false; // if true, the cellZone heats up because of the solar rad
        }
    }

SourceFiles
    absorptionEmissionSource.C

\*---------------------------------------------------------------------------*/
#ifndef absorptionEmissionSource_H
#define absorptionEmissionSource_H

#include "cellSetOption/cellSetOption.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
              Class absorptionEmissionSource Declaration
\*---------------------------------------------------------------------------*/

class absorptionEmissionSource
:
    public cellSetOption
{

protected:

        //- Name of the T equation field
        word TName_;

        //- Convert solar load into heat
        bool solarHeat_;

        //- Emission coeff. [1/m]
        scalar e_;

        //- Absorption coeff. [1/m]
        scalar a_;

        // Emission contribution [kg/(m s^3)]
        autoPtr<Function1<scalar>> E_;


    // Private Member Functions

        //- Calculate Ru for energy eqn.
        void calculateRu(scalarField& Ru);

        //- Calculate Rp for energy eqn.
        void calculateRp(scalarField& Rp);


public:

    //- Runtime type information
    TypeName("absorptionEmissionSource");


    // Constructors

        // - Construct from components
        absorptionEmissionSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~absorptionEmissionSource()
    {}

    // Member Functions

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames);

        virtual void addSup(fvMatrix<scalar>& eqn, const label fieldi);

        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Temperature equation sources
        void addSupT(fvMatrix<scalar>& eqn);

        //- Enthalpy equation sources
        virtual void addSupHe(fvMatrix<scalar>& eqn);

        //- ILambda equations absorption and emission terms
        virtual void addSupILambda(fvMatrix<scalar>& eqn);

        //- Isolar source
        void addSupIsolar(fvMatrix<scalar>& eqn);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const absorptionEmissionSource&) = delete;
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
