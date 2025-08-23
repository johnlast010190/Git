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
    (c) 2012-2013 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::fv::temperatureLimitsConstraint

Description
    Constraint for temperature to apply limits between minimum and maximum
    values

    Constraint described by:

        temperatureLimitsConstraintCoeffs
        {
            minimum     200;
            maximum     500;
        }


SourceFiles
    fvOption.C

\*---------------------------------------------------------------------------*/

#ifndef temperatureLimitsConstraint_H
#define temperatureLimitsConstraint_H

#include "cellSetOption/cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class temperatureLimitsConstraint Declaration
\*---------------------------------------------------------------------------*/

class temperatureLimitsConstraint
:
    public cellSetOption
{
protected:

    // Protected data

        //- Minimum temperature limit [K]
        scalar Tmin_;

        //- Maximum temperature limit [K]
        scalar Tmax_;


public:

    //- Runtime type information
    TypeName("temperatureLimitsConstraint");


    // Constructors

        //- Construct from components
        temperatureLimitsConstraint
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        temperatureLimitsConstraint
        (
            const temperatureLimitsConstraint&
        ) = delete;


    //- Destructor
    virtual ~temperatureLimitsConstraint()
    {}


    // Member Functions

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames) override;

        virtual bool alwaysApply() const;


        // Evaluate

            //- Correct the energy field
            virtual void correct(volScalarField& he) override;


        // I-O

            //- Write data
            virtual void writeData(Ostream&) const override;

            //- Read dictionary
            virtual bool read(const dictionary& dict) override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const temperatureLimitsConstraint&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
