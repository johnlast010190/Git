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
    (c) 2021-2024 Engys Ltd.
    (c) 2017 OpenCFD Ltd.
    (c) 2011-2023 OpenFOAM Foundation

Class
    Foam::Reaction

Description
    Simple extension of ThermoType to handle reaction kinetics in addition
    to the equilibrium thermodynamics already handled.

SourceFiles
    ReactionI.H
    Reaction.C

\*---------------------------------------------------------------------------*/

#ifndef Reaction_H
#define Reaction_H

#include "containers/HashTables/HashPtrTable/HashPtrTable.H"
#include "fields/Fields/scalarField/scalarField.H"
#include "db/typeInfo/typeInfo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "include/dummyThermo.H"
#include "materialModels/materialTables/materialTables.H"
#include "general/referenceFields/referenceFields.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"
#include "mixtures/coefficientMixture/coefficientMixture.H"
#include "matrices/simpleMatrix/simpleMatrix.H"
#include "primitives/Tuple2/Tuple2.H"
#include "reaction/reaction/reaction.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class ThermoType>
class Reaction;

template<class ThermoType>
inline Ostream& operator<<(Ostream&, const Reaction<ThermoType>&);


/*---------------------------------------------------------------------------*\
                           Class Reaction Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class Reaction
:
    public reaction,
    public ThermoType::thermoType
{
public:

    // Static data

        //- Default temperature limits of applicability of reaction rates
        static scalar TlowDefault, ThighDefault;


private:

    // Private data

        //- Temperature limits of applicability of reaction rates
        scalar Tlow_, Thigh_;

        //- Material tables pointer
        autoPtr<materialTables> mat_;

        //- Is it using new material library?
        mutable bool isMaterials_;

        //- Are the material model initialized?
        mutable bool isMatInit_;

        //- Mollar weights for the chemical equation
        mutable scalarField eqnWs_;

        //- Difference btween left/right side of moles in chemical equation
        mutable label nMoles_;

        //- Pointer to the gibs energy models
        mutable matScalarUTable gModels_;

        //- Pointer to the entropy models
        mutable matScalarUTable sModels_;

        //- Stoich coeffs
        mutable scalarField stoichCoeffs_;

        //- Phase name
        mutable word phaseName_;

        //- Chemistry coefficient mixture models
        mutable autoPtr<coefficientMixture> coeffMixture_;


    // Private Member Functions

        //- Construct reaction thermo
        void setThermo(const PtrList<ThermoType>& speciesThermo);

        //- Disallow default bitwise assignment
        void operator=(const Reaction<ThermoType>&) = delete;

        //- Initialize material models
        void initMaterialModels() const;


public:

    //- Runtime type information
    TypeName("Reaction");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            Reaction,
            dictionary,
            (
                const speciesTable& species,
                const PtrList<ThermoType>& speciesThermo,
                const dictionary& dict
            ),
            (species, speciesThermo, dict)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            Reaction,
            objectRegistry,
            (
                const speciesTable& species,
                const PtrList<ThermoType>& speciesThermo,
                const objectRegistry& obr,
                const dictionary& dict
            ),
            (species, speciesThermo, obr, dict)
        );


    // Constructors

        //- Construct from components
        Reaction
        (
            const speciesTable& species,
            const PtrList<ThermoType>& speciesThermo,
            const List<specieCoeffs>& lhs,
            const List<specieCoeffs>& rhs
        );

        //- Construct as copy given new speciesTable
        Reaction
        (
            const Reaction<ThermoType>&,
            const speciesTable& species
        );

        //- Construct from dictionary
        // NOTE: Initialise thermo is used by solidReaction where there is no
        // need of setting a lhs - rhs thermo type for each reaction. This is
        // needed for mechanism with reversible reactions
        Reaction
        (
            const speciesTable& species,
            const PtrList<ThermoType>& speciesThermo,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual autoPtr<Reaction<ThermoType>> clone() const = 0;

        //- Construct and return a clone with new speciesTable
        virtual autoPtr<Reaction<ThermoType>> clone
        (
            const speciesTable& species
        ) const = 0;


    // Selectors

        //- Return a pointer to new patchField created on freestore
        //  from dictionary
        static autoPtr<Reaction<ThermoType>> New
        (
            const speciesTable& species,
            const PtrList<ThermoType>& speciesThermo,
            const dictionary& dict
        );

        //- Return a pointer to new patchField created on freestore
        //  from objectRegistry and dictionary
        static autoPtr<Reaction<ThermoType>> New
        (
            const speciesTable& species,
            const PtrList<ThermoType>& speciesThermo,
            const objectRegistry& obr,
            const dictionary& dict
        );


    //- Destructor
    virtual ~Reaction()
    {}


    // Member Functions

        // Access

            //- Return the name of the reaction
            using reaction::name;

            //- Return the lower temperature limit for the reaction
            inline scalar Tlow() const;

            //- Return the upper temperature limit for the reaction
            inline scalar Thigh() const;

        // Hooks

            //- Pre-evaluation hook
            virtual void preEvaluate() const = 0;

            //- Post-evaluation hook
            virtual void postEvaluate() const = 0;


        // Reaction rate coefficients

            //- Concentration powers
            void C
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalar& Cf,
                scalar& Cr
            ) const;

            //- Net reaction rate
            scalar omega
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalar& omegaf,
                scalar& omegar
            ) const;

            //- The net reaction rate for each species involved
            void dNdtByV
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalarField& dNdtByV,
                const bool reduced,
                const List<label>& c2s,
                const label Nsi0
            ) const;


        // Reaction rate coefficients

            //- Forward rate constant
            virtual scalar kf
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const = 0;

            //- Reverse rate constant from the given forward rate constant
            virtual scalar kr
            (
                const scalar kf,
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const = 0;

            //- Reverse rate constant.
            //  Note this evaluates the forward rate constant and divides by the
            //  equilibrium constant
            virtual scalar kr
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const = 0;


        // Jacobian coefficients

            //- Temperature derivative of forward rate
            virtual scalar dkfdT
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const = 0;

            //- Temperature derivative of reverse rate
            virtual scalar dkrdT
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                const scalar dkfdT,
                const scalar kr
            ) const = 0;

            //- Does this reaction have concentration-dependent rate constants?
            virtual bool hasDkdc() const = 0;

            //- Concentration derivative of forward rate
            virtual void dkfdc
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalarField& dkfdc
            ) const = 0;

            //- Concentration derivative of reverse rate
            virtual void dkrdc
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                const scalarField& dkfdc,
                const scalar kr,
                scalarField& dkrdc
            ) const = 0;

            //- Derivative of the net reaction rate for each species involved
            //  w.r.t. the concentration and temperature
            void ddNdtByVdcTp
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalarField& dNdtByV,
                scalarSquareMatrix& ddNdtByVdcTp,
                const bool reduced,
                const List<label>& c2s,
                const label csi0,
                const label Tsi,
                scalarField& cTpWork0,
                scalarField& cTpWork1
            ) const;

            //- Equilibrium constant i.t.o. molar concentration
            //  = PIi(ci/cStd)^nui
            //  For low pressures (where the gas mixture is near perfect)
            //  Kc = Kp(Pstd/(RR*T))^nu
            scalar Kc(const scalar p, const scalar T) const;

            //- Derivative of B (according to Niemeyer et al.)
            //  w.r.t. temperature
            scalar dKcdTbyKc(const scalar p, const scalar T) const;

        //- Write
        virtual void write(Ostream&) const;


    // Ostream Operator

        friend Ostream& operator<< <ThermoType>
        (
            Ostream&,
            const Reaction<ThermoType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "reaction/Reactions/Reaction/ReactionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "reaction/Reactions/Reaction/Reaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
