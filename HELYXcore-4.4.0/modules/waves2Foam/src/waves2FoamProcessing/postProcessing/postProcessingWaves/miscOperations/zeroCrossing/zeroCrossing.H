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
    (c) held by original author
    (c) 2024 Engys Ltd.

Class
    Foam::zeroCrossing

Description
    Perform an ensemble average of a given data stream.

SourceFiles
    zeroCrossing.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef zeroCrossing_H
#define zeroCrossing_H

#include "postProcessing/postProcessingWaves/postProcessingWaves.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class zeroCrossing
:
    public postProcessingWaves
{
protected:

    // Protected data

        labelList indices_;
        scalar deltaT_;

        word crossingType_;

        IOdictionary dataDict_;


private:

    // Private Member Functions

        void evaluateScalar();

        void writeScalar
        (
            const List<scalarField>&,
            const List<scalarField>&,
            const List<scalarField>&,
            const List<scalarField>&,
            const word
        ) const;

        void crossingAnalysis
        (
            const scalarField&,
            const List<scalarField>&,
            List<scalarField>&,
            List<scalarField>&,
            List<scalarField>&,
            List<scalarField>&,
            const label
        ) const;


public:

    //- Runtime type information
    TypeName("zeroCrossing");


    // Constructors

        //- Construct from components
        zeroCrossing
        (
            const Time&,
            const dictionary&,
            const word&
        );

        //- Disallow default bitwise copy construction
        zeroCrossing(const zeroCrossing&) = delete;


    //- Destructor
    virtual ~zeroCrossing();


    // Member Functions

        virtual void evaluate();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const zeroCrossing&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
