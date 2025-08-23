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
    (c) 2022-2023 OpenFOAM Foundation
    (c) 2024 Engys Ltd

Class
    Foam::fv::phaseTurbulenceStabilisation

Description
    Phase turbulence stabilisation

    In the limit of a phase-fraction->0 the turbulence properties cannot be
    obtained from the phase turbulence model, coupling to the other phase/phases
    is required.  The phaseTurbulenceStabilisation fvModel stabilises the phase
    turbulence properties by adding transfer terms from the corresponding
    properties of the other phases when the phase-fraction is less than the
    specified \c alphaInversion.  This implementation is a generalisation of
    the approach used in the Foam::RASModels::LaheyKEpsilon and
    Foam::RASModels::continuousGasKEpsilon models to handle phase-inversion and
    free-surface flow and can be used with any combination of RAS turbulence
    models.

    To stabilise the solution of the phase turbulence equations \c
    alphaInversion can be set to a small value e.g. 1e-2, but unless the phase
    turbulence model is specifically designed to handle phase-inversion and both
    continuous and dispersed regimes it may be useful to set \c alphaInversion
    to a higher value, corresponding to the phase-fraction at which transition
    from continuous to dispersed happens and effectively use the turbulence
    properties of the other phase when the phase is dispersed.  This is of
    course an approximation to the real system and if accurate handling of both
    the continuous and dispersed phase regimes is required specially developed
    models should be used.

Usage
    Example usage:
    \verbatim
    phaseTurbulenceStabilisation
    {
        type    phaseTurbulenceStabilisation;
        phase   air;
        alphaInversion  0.1;
    }
    \endverbatim

SourceFiles
    phaseTurbulenceStabilisation.C

\*---------------------------------------------------------------------------*/

#ifndef phaseTurbulenceStabilisation_H
#define phaseTurbulenceStabilisation_H

#include "cfdTools/general/fvOptions/fvOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class eulerianPhaseSystem;
class turbulenceModel;

namespace fv
{

/*---------------------------------------------------------------------------*\
                Class phaseTurbulenceStabilisation Declaration
\*---------------------------------------------------------------------------*/

class phaseTurbulenceStabilisation
:
    public option
{
    // Private Data

        //- The name of the phase
        word phaseName_;

        //- Phase-fraction below which turbulence property blending is applied
        dimensionedScalar alphaInversion_;

        //- Reference to the phase
        const eulerianPhaseSystem* phaseSystem_;

        //- Reference to the mixture turbulence model
        const turbulenceModel* turbulence_;


    // Private Member Functions

        //- Add contribution to phase psi equation
        void addAlphaRhoSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const volScalarField& field,
            fvMatrix<scalar>& eqn,
            tmp<volScalarField>
            (turbulenceModel::*psi)() const
        ) const;


public:

    //- Runtime type information
    TypeName("phaseTurbulenceStabilisation");


    // Constructors

        //- Construct from explicit source name and mesh
        phaseTurbulenceStabilisation
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        phaseTurbulenceStabilisation
        (
            const phaseTurbulenceStabilisation&
        ) = delete;


    // Member Functions

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames) override;


        // Sources

            //- Add contribution to phase k, epsilon or omega equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            ) override;


        // Mesh changes

            //- Update for mesh motion
            virtual bool movePoints()
            {
                return true;
            }

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&){}

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&){}

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&){}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const phaseTurbulenceStabilisation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
