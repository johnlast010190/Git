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
    (c) 2023-2025 Engys Ltd.

Class
    Foam::multiphaseThermo

SourceFiles
    multiphaseThermo.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseThermo_H
#define multiphaseThermo_H

#include "fields/volFields/volFields.H"
#include "rhoThermo/rhoThermo.H"
#include "materialModels/materialTables/materialTables.H"
#include "mixtures/phaseVolumeFractions/phaseVolumeFractions.H"
#include "phasePairKey/phasePairKey.H"
#include "binaryPhaseModels/surfaceTensionModels/surfaceTensionModel/surfaceTensionModel.H"
#include "binaryPhaseModels/phaseChangeModels/phaseChangeModel/phaseChangeModel.H"
#include "solidThermo/solidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class multiphaseThermo Declaration
\*---------------------------------------------------------------------------*/

class multiphaseThermo
:
    public rhoThermo::composite
{
public:

    // Public typedefs

        typedef
            HashTable
            <
                autoPtr<surfaceTensionModel>,
                phasePairKey,
                phasePairKey::hash
            >
            surfaceTensionModelTable;

        typedef
            HashTable
            <
                autoPtr<phaseChangeModel>,
                phasePairKey,
                phasePairKey::hash
            >
            phaseChangeModelTable;


private:

    // Private Data

        //- Mesh reference.
        const fvMesh& mesh_;

        //- Reference to material models
        materialTables& materials_;

        //- Pointer list for phases thermos
        PtrList<basicThermo> thermos_;

        //- Phase list
        const wordList phases_;

        //- Pointer to volume fraction class
        phaseVolumeFractions* fractions_;

        //- Surface tension models
        surfaceTensionModelTable surfaceTensionModels_;

        //- Hash table for sigma interphase pairs
        typedef
            HashPtrTable<Function1<scalar>, PairKey<word>, PairKey<word>::hash>
            sigmaTable;

        //- Phase-change models
        phaseChangeModelTable phaseChangeModels_;

        //- Helper class to create Function1's from stream for HashPtrTable
        class iNewSigma
        {
        public:

            iNewSigma()
            {}

            autoPtr<Function1<scalar>> operator()
            (
                const PairKey<word>& key, Istream& is
            ) const
            {
                word keyName = "sigma";
                token firstToken(is);
                if (firstToken.isScalar())
                {
                    dictionary d;
                    d.add(keyName, firstToken);
                    return
                        autoPtr<Function1<scalar>>
                        (
                            Function1<scalar>::New(keyName, d)
                        );
                }
                is.putBack(firstToken);
                autoPtr<entry> e = entry::New(is);
                dictionary d;
                d.add(e);
                return
                    autoPtr<Function1<scalar>>
                    (
                        Function1<scalar>::New(keyName, d)
                    );
            }
        };

        //- Legacy - surface tension
        sigmaTable sigmas_;


    // Private Member Functions

        //- Construct material tables
        materialTables& matLookupOrConstruct
        (
            const objectRegistry& obr,
            const dictionary& dict
        );

        //- Return boundary field for different combination of
        // temperature/pressure
        tmp<scalarField> tmpBoundaryField
        (
            const scalarField& p,
            const scalarField& T,
            const label patchi,
            const word& modelName
        ) const;


public:

    //- Runtime type information
    TypeName("multiphase");


    // Constructors

        //- Construct from object registry
        // (phase name not applicable but necessary
        // to add to run time selection tables)
        multiphaseThermo
        (
            const objectRegistry& obr,
            const word& phaseName = word::null
        );


    // Selectors

        //- Generic lookup or create for any thermo type
        //  Returns the sub-thermo for the specified phaseName if given
        static basicThermo& lookupOrCreate
        (
            const objectRegistry&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~multiphaseThermo()
    {}


    // Member Functions

        //- Additional initialisation - called by parent BasicThermo after
        //  construction
        virtual void init() override;

        //- Update properties
        virtual void correct() override;

        //- Update densities for given pressure change
        void correctRho(const volScalarField& dp) override;

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        virtual bool incompressible() const override;

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        virtual bool isochoric() const override;

        //- Thermos for phases
        const PtrList<basicThermo>& thermos() const
        {
            return thermos_;
        }

        //- Reference to thermos for phases
        PtrList<basicThermo>& thermos()
        {
            return thermos_;
        }

        //- Thermos for phase
        const basicThermo& thermos(const word& phaseName) const
        {
            const label idx = phases_.find(phaseName);
            if (idx < 0)
            {
                FatalErrorInFunction
                    << "Phase " + phaseName + " not found. Available phases: "
                    << phases_ << nl << exit(FatalError);
            }
            return thermos_[idx];
        }

        //- Return solid thermo in multiphase
        const solidThermo& sThermo() const
        {
            label nSolidThermos = 0;
            const solidThermo* lastSolidThermo = nullptr;
            forAll(thermos_, i)
            {
                const basicThermo* thermo = &thermos_[i];
                if (isA<solidThermo>(*thermo))
                {
                    ++nSolidThermos;
                    lastSolidThermo = dynamic_cast<const solidThermo*>(thermo);
                }
            }
            if (nSolidThermos != 1)
            {
                FatalErrorInFunction
                    << "Only one solidThermo is allowed" << abort(FatalError);
                NotImplemented;
            }
            return *lastSolidThermo;
        }

        //- Thermos for phases
        const PtrList<volScalarField>& alphas() const
        {
            return fractions_->alphas();
        }

        //- Thermos for phases
        PtrList<volScalarField>& alphas()
        {
            return fractions_->alphas();
        }

        //- Thermos for pahses
        PtrList<volScalarField>& dgdts() const
        {
            return const_cast<PtrList<volScalarField>&>(fractions_->dgdts());
        }

        const wordList& phases() const
        {
            return phases_;
        }

        phaseVolumeFractions& fractions()
        {
            HELYX_ASSERT(fractions_)
            {
                FatalErrorInFunction
                    << "Volume fractions is uninitialised"
                    << exit(FatalError);
            }
            return *fractions_;
        }

        // Derived thermodynamic properties

            //- Return the phase-averaged reciprocal Cv
            tmp<volScalarField> rCv() const;


        // Phase-pair data

            //- Gather and create binary models
            template<class modelType>
            void createBinaryModels
            (
                HashTable
                <
                    autoPtr<modelType>,
                    phasePairKey,
                    phasePairKey::hash
                >& models,
                const bool ordered,
                const word& orderedKey = "in"
            );

            //- Return the subdict containing the stipulated keyword
            //  This is to allow multiple intersecting dictionaries with wild-
            //  cards to co-exist where different models are specified in each
            //  If modelName == word::null, don't check for the stipulated
            //  keyword.
            //  If doRegex == false, don't do wildcard matching
            const dictionary* phasePairDictContaining
            (
                const word& modelName,
                const word& phase1,
                const word& phase2,
                const bool ordered,
                const bool doRegex=true
            );

            //- Same as above but with dictionary given
            static dictionary* phasePairDictContaining
            (
                dictionary& dict,
                const word& modelName,
                const word& phase1,
                const word& phase2,
                const bool ordered,
                const bool doRegex=true
            );

            surfaceTensionModelTable& surfaceTensionModels()
            {
                return surfaceTensionModels_;
            }

            phaseChangeModelTable& phaseChangeModels()
            {
                return phaseChangeModels_;
            }

            const surfaceTensionModel* surfaceTensionModels
            (
                const word& phase1, const word& phase2
            ) const;

            const phaseChangeModel* phaseChangeModels
            (
                const word& phase1, const word& phase2
            ) const;


    // TODO Functions and data to be moved

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

        //- Dot product of face unit interface normal flux
        // and Face area vectors
        tmp<surfaceVectorField> nHatfv
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        //- Face unit interface normal flux
        tmp<surfaceScalarField> nHatf
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        //- Function correcting contact angle
        void correctContactAngle
        (
            const volVectorField& U,
            const volScalarField& alpha1,
            const volScalarField& alpha2,
            surfaceVectorField::Boundary& nHatb
        ) const;

        //- Simple expression for curvature
        tmp<volScalarField> K
        (
            const volVectorField& U,
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        //- Surface tension force between two phases
        tmp<surfaceScalarField> surfaceTensionForce
        (
            const volVectorField& U,
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        //- Surface tension force acting on a particular phase
        tmp<surfaceScalarField> surfaceTensionForce
        (
            const volVectorField& U,
            const volScalarField& alpha1
        ) const;

        //- Total surface tension force
        tmp<surfaceScalarField> surfaceTensionForce
        (
            const volVectorField& U
        ) const;

        //- Indicator of the proximity of the interface
        //  Field values are 1 near and 0 away for the interface.
        tmp<volScalarField> nearInterface() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multiphaseThermoTemplates.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
