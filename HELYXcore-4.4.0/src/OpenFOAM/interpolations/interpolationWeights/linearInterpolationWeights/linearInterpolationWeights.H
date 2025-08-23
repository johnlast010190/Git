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
    (c) 2012-2019 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::linearInterpolationWeights

Description

SourceFiles
    linearInterpolationWeights.C

\*---------------------------------------------------------------------------*/

#ifndef linearInterpolationWeights_H
#define linearInterpolationWeights_H

#include "interpolations/interpolationWeights/interpolationWeights/interpolationWeights.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class linearInterpolationWeights Declaration
\*---------------------------------------------------------------------------*/

class linearInterpolationWeights
:
    public interpolationWeights
{

private:

    // Private data

        //- Cached index in samples from previous invocation
        mutable label index_;


public:

    //- Runtime type information
    TypeName("linear");

    // Constructors

        //- Construct from components
        linearInterpolationWeights
        (
            const scalarField& samples
        );


    //- Destructor
    virtual ~linearInterpolationWeights()
    {}


    // Member Functions

        //- Calculate weights and indices to calculate t from samples.
        //  Returns true if indices changed.
        virtual bool valueWeights
        (
            const scalar t,
            labelList& indices,
            scalarField& weights
        ) const override;

        //- Calculate weights and indices to calculate integrand of t1..t2
        //  from samples. Returns true if indices changed.
        virtual bool integrationWeights
        (
            scalar t1,
            scalar t2,
            labelList& indices,
            scalarField& weights
        ) const override;

        //- Calculate weights and indices to calculate derivative of t1..t2
        //  from samples. Returns true if indices changed.
        virtual bool derivationWeights
        (
            const scalar x,
            labelList& indices,
            scalarField& weights
        ) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
