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
    (c) 2020-2024 Engys Ltd.

Class
    Foam::fv::settlingVelocityTakacs

Description

    solverObject to create, store and compute a settling flux to be used
    in an additional passive scalar transport solverObject.

SourceFiles
    settlingVelocityTakacs.C

\*---------------------------------------------------------------------------*/
#ifndef settlingVelocityTakacs_H
#define settlingVelocityTakacs_H

#include "cellSetOption/cellSetOption.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class settlingVelocityTakacs Declaration
\*---------------------------------------------------------------------------*/

class settlingVelocityTakacs
:
    public cellSetOption
{
    // Private data

        //- Concentration field name
        word XName_;

        //- Flux field name
        word phiName_;

        //- Phase field name
        word phaseName_;

        //- Model parameter
        dimensionedVector V0_;
        dimensionedVector V0t_;
        dimensionedScalar rh_;
        dimensionedScalar rp_;
        dimensionedScalar rt_;
        dimensionedScalar Xns_;
        dimensionedScalar X34_;

        //- Scaling parameter to get X in [g/L]
        scalar fac_;
        scalar blend_;

        //- use constant drift y/n
        bool constant_;

        //- limit phasic drift
        bool limitBulk_;
        scalar threshold_;
        scalar limitHeight_;


    // Private Member Functions

        //- compute settling flux
        tmp<surfaceScalarField> computeFlux(const fvMesh& mesh);

        //- limit to selected cells
        void limitCellSet(volVectorField& Udrift);


public:

    //- Runtime type information
    TypeName("settlingVelocityTakacs");


    // Constructors

        //- Construct from components
        settlingVelocityTakacs
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        settlingVelocityTakacs(const settlingVelocityTakacs&) = delete;


    // Member Functions

        //- Initialise after all fvOptions created
        virtual bool initialise() override;

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames) override;


        // Evaluate

            //- Add settling advection term to scalar transport equation
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldi
            ) override;

            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            ) override;

        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict) override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const settlingVelocityTakacs&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
