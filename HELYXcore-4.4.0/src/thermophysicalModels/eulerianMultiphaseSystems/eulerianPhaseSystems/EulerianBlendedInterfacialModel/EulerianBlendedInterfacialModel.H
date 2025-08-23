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
    (c) 2014-2016 OpenFOAM Foundation

Class
    Foam::EulerianBlendedInterfacialModel

Description

SourceFiles
    EulerianBlendedInterfacialModel.C

\*---------------------------------------------------------------------------*/

#ifndef BlendedInterfacialModel_H
#define BlendedInterfacialModel_H

#include "eulerianBlendingMethods/eulerianBlendingMethod/eulerianBlendingMethod.H"
#include "../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "../eulerianPhasePair/orderedEulerianPhasePair/orderedEulerianPhasePair.H"

#include "fields/GeometricFields/geometricZeroField/geometricZeroField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class EulerianBlendedInterfacialModel Declaration
\*---------------------------------------------------------------------------*/

template<class ModelType>
class EulerianBlendedInterfacialModel
{
    // Private data

        //- Reference to phase 1
        const eulerianPhaseModel& phase1_;

        //- Reference to phase 2
        const eulerianPhaseModel& phase2_;

        //- Blending model
        const eulerianBlendingMethod& blending_;

        //- Model for region with no obvious dispersed phase
        autoPtr<ModelType> model_;

        //- Model for dispersed phase 1 in continuous phase 2
        autoPtr<ModelType> model1In2_;

        //- Model for dispersed phase 2 in continuous phase 1
        autoPtr<ModelType> model2In1_;

        //- If true set coefficients and forces to 0 at fixed-flux BCs
        bool correctFixedFluxBCs_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        EulerianBlendedInterfacialModel(const EulerianBlendedInterfacialModel<ModelType>&);

        //- Disallow default bitwise assignment
        void operator=(const EulerianBlendedInterfacialModel<ModelType>&);

        //- Correct coeff/value on fixed flux boundary conditions
        template<class GeometricField>
        void correctFixedFluxBCs(GeometricField& field) const;


public:

    // Constructors

        //- Construct from two phases, blending method and three models
        EulerianBlendedInterfacialModel
        (
            const eulerianPhaseModel& phase1,
            const eulerianPhaseModel& phase2,
            const eulerianBlendingMethod& blending,
            autoPtr<ModelType> model,
            autoPtr<ModelType> model1In2,
            autoPtr<ModelType> model2In1,
            const bool correctFixedFluxBCs = true
        );


        //- Construct from the model table, dictionary and pairs
        EulerianBlendedInterfacialModel
        (
            const eulerianPhasePair::dictTable& modelTable,
            const eulerianBlendingMethod& blending,
            const eulerianPhasePair& pair,
            const orderedEulerianPhasePair& pair1In2,
            const orderedEulerianPhasePair& pair2In1,
            const bool correctFixedFluxBCs = true
        );


    //- Destructor
    ~EulerianBlendedInterfacialModel();


    // Member Functions

        //- Return true if a model is specified for the supplied phase
        bool hasModel(const eulerianPhaseModel& phase) const;

        //- Return the model for the supplied phase
        const ModelType& model(const eulerianPhaseModel& phase) const;

        //- Return the sign of the explicit value for the supplied phase
        scalar sign(const eulerianPhaseModel& phase) const;

        //- Return the blended force coefficient
        tmp<volScalarField> K() const;

        //- Return the blended force coefficient with a specified residual alpha
        tmp<volScalarField> K(const scalar residualAlpha) const;

        //- Return the face blended force coefficient
        tmp<surfaceScalarField> Kf() const;

        //- Return the blended force
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>> F() const;

        //- Return the face blended force
        tmp<surfaceScalarField> Ff() const;

        //- Return the blended diffusivity
        tmp<volScalarField> D() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "EulerianBlendedInterfacialModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
