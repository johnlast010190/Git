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
    (c) 2024-2025 Engys Ltd

Class
    Foam::phaseChangeModel

Description

SourceFiles
    phaseChangeModel.C
    newPhaseChangeModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseChangeModel2_H
#define phaseChangeModel2_H

#include "db/typeInfo/typeInfo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "fields/volFields/volFields.H"
#include "dimensionedTypes/dimensionedScalar/dimensionedScalar.H"
#include "memory/autoPtr/autoPtr.H"
#include "primitives/Pair/Pair.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "phasePairKey/phasePairKey.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class multiphaseThermo;

/*---------------------------------------------------------------------------*\
                     Class phaseChangeModel Declaration
\*---------------------------------------------------------------------------*/

class phaseChangeModel
{
protected:

    // Protected data

        const fvMesh& mesh_;

        //- Phase pair
        const phasePairKey& pair_;

        const multiphaseThermo& thermo_;

        //- Indices of the two phases in list of thermos
        const label phaseIdx1_;
        const label phaseIdx2_;

        //- Saturation vapour pressure
        dimensionedScalar pSat_;

        //- Saturation vapour pressure table values
        //  (to gradually activate cavitation)
        autoPtr<Function1<scalar>> pSatTable_;


public:

    //- Runtime type information
    TypeName("phaseChangeModel");


    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        phaseChangeModel,
        multiphase,
        (
            const dictionary& dict,
            const phasePairKey& pair,
            const fvMesh& mesh
        ),
        (dict, pair, mesh)
    );

    // Constructors

        //- Construct from a dictionary and a phase pair
        phaseChangeModel
        (
            const dictionary& dict,
            const phasePairKey& pair,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        phaseChangeModel(const phaseChangeModel&) = delete;


    //- Destructor
    virtual ~phaseChangeModel()
    {}


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<phaseChangeModel> New
        (
            const dictionary& dict,
            const phasePairKey& pair,
            const fvMesh& mesh
        );


    // Member Functions

        const volScalarField& alpha1() const;

        const volScalarField& alpha2() const;

        tmp<volScalarField> rho1() const;

        tmp<volScalarField> rho2() const;

        //- Return const-access to the saturation vapour pressure
        const dimensionedScalar& pSat() const
        {
            return pSat_;
        }

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        virtual Pair<tmp<volScalarField>> mDotAlphal() const = 0;

        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply (p - pSat)
        virtual Pair<tmp<volScalarField>> mDotP() const = 0;

        //- Return the volumetric condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        Pair<tmp<volScalarField>> vDotAlphal() const;

        //- Return the volumetric condensation and vaporisation rates as
        //  coefficients to multiply (p - pSat)
        Pair<tmp<volScalarField>> vDotP() const;

        //- Correct the phaseChange model
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const phaseChangeModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
