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
    (c) 2014-2015 OpenFOAM Foundation

Class
    Foam::eulerianPhasePair

Description

SourceFiles
    eulerianPhasePair.C

\*---------------------------------------------------------------------------*/

#ifndef phasePair_H
#define phasePair_H

#include "../../eulerianPhaseModel/eulerianPhaseModel/eulerianPhaseModel.H"
#include "phasePairKey/phasePairKey.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class eulerianPhasePair Declaration
\*---------------------------------------------------------------------------*/

class eulerianPhasePair
:
    public phasePairKey
{
public:

    // Hash table types

        //- Dictionary hash table
        typedef HashTable<dictionary, phasePairKey, phasePairKey::hash>
            dictTable;

        //- Scalar hash table
        typedef HashTable<scalar, phasePairKey, phasePairKey::hash>
            scalarTable;


private:

    // Private data

        //- Phase 1
        const eulerianPhaseModel& phase1_;

        //- Phase 2
        const eulerianPhaseModel& phase2_;

        //- Gravitational acceleration if found or nullptr otherwise
        const uniformDimensionedVectorField* gPtr_;


    // Private member functions

        // Etvos number for given diameter
        tmp<volScalarField> EoH(const volScalarField& d) const;


public:

    // Constructors

        //- Construct from two phases and gravity
        eulerianPhasePair
        (
            const eulerianPhaseModel& phase1,
            const eulerianPhaseModel& phase2,
            const bool ordered = false
        );


    //- Destructor
    virtual ~eulerianPhasePair();


    // Member Functions

        //- Dispersed phase
        virtual const eulerianPhaseModel& dispersed() const;

        //- Continuous phase
        virtual const eulerianPhaseModel& continuous() const;

        //- Average density
        tmp<volScalarField> rho() const;

        //- Relative velocity magnitude
        tmp<volScalarField> magUr() const;

        //- Relative velocity
        tmp<volVectorField> Ur() const;

        //- Reynolds number
        tmp<volScalarField> Re() const;

        //- Prandtl number
        tmp<volScalarField> Pr() const;

        //- Eotvos number
        tmp<volScalarField> Eo() const;

        //- Eotvos number based on hydraulic diameter type 1
        tmp<volScalarField> EoH1() const;

        //- Eotvos number based on hydraulic diameter type 2
        tmp<volScalarField> EoH2() const;

        //- Surface tension coefficient
        tmp<volScalarField> sigma() const;

        //- Morton Number
        tmp<volScalarField> Mo() const;

        //- Takahashi Number
        tmp<volScalarField> Ta() const;

        //- Aspect ratio
        virtual tmp<volScalarField> E() const;

        // Access

            // Phase 1
            inline const eulerianPhaseModel& phase1() const;

            // Phase 2
            inline const eulerianPhaseModel& phase2() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "eulerianPhasePairI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
