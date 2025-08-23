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
    (c) 2015-2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd

Class
    Foam::eulerianPhaseSystem

Description
    Class to represent a system of phases and model interfacial transfers
    between them.

SourceFiles
    eulerianPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef phaseSystem_H
#define phaseSystem_H

#include "db/regIOobject/regIOobject.H"

#include "../eulerianPhaseModel/eulerianPhaseModel/eulerianPhaseModel.H"
#include "../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "../eulerianPhasePair/orderedEulerianPhasePair/orderedEulerianPhasePair.H"
#include "containers/HashTables/HashPtrTable/HashPtrTable.H"
#include "containers/Dictionaries/PtrListDictionary/PtrListDictionary.H"

#include "cfdTools/general/fvOptions/fvOptions.H"

#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fvMatrices/fvMatricesFwd.H"

#include "multiphaseThermo/multiphaseThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class eulerianBlendingMethod;
template<class modelType> class EulerianBlendedInterfacialModel;
class eulerianAspectRatioModel;

/*---------------------------------------------------------------------------*\
                         Class eulerianPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

class eulerianPhaseSystem
:
    public regIOobject
{
public:

    // Public typedefs

        typedef
            HashPtrTable
            <
                volScalarField,
                phasePairKey,
                phasePairKey::hash
            >
            KdTable;

        typedef
            HashPtrTable
            <
                volScalarField,
                phasePairKey,
                phasePairKey::hash
            >
            VmTable;

        typedef
            HashPtrTable
            <
                fvVectorMatrix,
                word,
                string::hash
            >
            momentumTransferTable;

        typedef
            HashPtrTable
            <
                fvScalarMatrix,
                word,
                string::hash
            >
            heatTransferTable;

        typedef
            HashPtrTable
            <
                fvScalarMatrix,
                word,
                string::hash
            >
            massTransferTable;

        typedef PtrListDictionary<eulerianPhaseModel> phaseModelList;


protected:

    // Protected typedefs

        typedef
            HashTable<dictionary, phasePairKey, phasePairKey::hash>
            dictTable;

        typedef
            HashTable<autoPtr<eulerianPhasePair>, phasePairKey, phasePairKey::hash>
            phasePairTable;

        typedef
            HashTable<autoPtr<eulerianBlendingMethod>, word, word::hash>
            blendingMethodTable;

        typedef
            HashTable
            <
                autoPtr<eulerianAspectRatioModel>,
                phasePairKey,
                phasePairKey::hash
            >
            aspectRatioModelTable;


    // Protected data

        //- Reference to the mesh
        const fvMesh& mesh_;

        multiphaseThermo& thermo_;

        //- Phase models
        phaseModelList phaseModels_;

        //- Phase pairs
        phasePairTable phasePairs_;

        //- Total volumetric flux
        surfaceScalarField phiv_;

        //- Rate of change of pressure
        volScalarField dpdt_;

        //- Blending methods
        blendingMethodTable blendingMethods_;


        // Sub Models

            //- Aspect ratio models
            aspectRatioModelTable aspectRatioModels_;


    // Protected member functions

        //- Calculate and return the mixture flux
        tmp<surfaceScalarField> calcPhiv
        (
            const phaseModelList& eulerianPhaseModels
        ) const;

        //- Generate pairs
        void generatePairs();

        //- Generate sub-model tables
        template<class modelType>
        void createSubModels
        (
            HashTable
            <
                autoPtr<modelType>,
                phasePairKey,
                phasePairKey::hash
            >& models,
            const bool ordered,
            const word& phaseName=word::null
        );

        //- Generate pairs and blended sub-model tables
        template<class modelType>
        void createSubModels
        (
            HashTable
            <
                autoPtr<EulerianBlendedInterfacialModel<modelType>>,
                phasePairKey,
                phasePairKey::hash
            >& models,
            const word& phaseName=word::null
        );

        //- Generate per-phase sub-model tables
        template<class modelType>
        void createSubModels
        (
            HashTable
            <
                HashTable<autoPtr<EulerianBlendedInterfacialModel<modelType>>>,
                phasePairKey,
                phasePairKey::hash
            >& models
        );


public:

    //- Runtime type information
    TypeName("phaseSystem");


    // Constructors

        //- Construct from fvMesh
        eulerianPhaseSystem(const fvMesh& mesh, const dictionary&);


    //- Destructor
    virtual ~eulerianPhaseSystem();


    // Member Functions

        //- Constant access the mesh
        inline const fvMesh& mesh() const;

        const multiphaseThermo& thermo() const
        {
            return thermo_;
        }

        //- Constant access the phase models
        inline const phaseModelList& phases() const;

        //- Access the phase models
        inline phaseModelList& phases();

        //- Constant access the phase pairs
        inline const phasePairTable& phasePairs() const;

        //- Return volume fractions from thermo model
        inline const PtrList<volScalarField>& alphas() const;

        //- Return volume fractions from thermo model
        inline PtrList<volScalarField>& alphas();

        //- Return the mixture density
        tmp<volScalarField> rho() const;

        //- Return the mixture velocity
        tmp<volVectorField> U() const;

        //- Return the mixture effective viscosity
        tmp<volScalarField> muEff() const;

        //- Constant access the mixture flux
        inline const surfaceScalarField& phiv() const;

        //- Access the mixture flux
        inline surfaceScalarField& phiv();

        //- Constant access the rate of change of the pressure
        inline const volScalarField& dpdt() const;

        //- Access the rate of change of the pressure
        inline volScalarField& dpdt();

        //- Return the aspect-ratio
        tmp<volScalarField> E(const phasePairKey& key) const;

        //- Optional FV-options
        inline fv::options& fvOptions() const;

        //- Access a sub model between a phase pair
        template<class modelType>
        const modelType& lookupSubModel(const eulerianPhasePair& key) const;

        //- Access a sub model between two phases
        template<class modelType>
        const modelType& lookupSubModel
        (
            const eulerianPhaseModel& dispersed,
            const eulerianPhaseModel& continuous
        ) const;

        //- Solve for the phase fractions
        virtual void solve();

        //- Correct the fluid properties other than the thermo and turbulence
        virtual void correct();

        //- Correct the kinematics
        virtual void correctKinematics();

        //- Correct the thermodynamics
        virtual void correctThermo();

        //- Correct the turbulence
        virtual void correctTurbulence();

        //- Correct the energy transport e.g. alphat
        virtual void correctEnergyTransport();

        virtual bool read();

        virtual bool writeData(Ostream&) const
        {
            return true;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "eulerianPhaseSystemI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "eulerianPhaseSystemTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
