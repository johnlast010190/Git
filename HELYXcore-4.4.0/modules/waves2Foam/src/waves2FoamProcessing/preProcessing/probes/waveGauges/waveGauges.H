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
    Foam::waveGauges

Description
    Writing wave gauges depending on a given point distribution.

SourceFiles
    waveGauges.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef waveGauges_H
#define waveGauges_H

#include "fvMesh/fvMesh.H"
#include "cfdTools/general/include/fvCFD.H"
#include "fields/volFields/volFields.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"

#include "preProcessing/probes/pointDistributions/pointDistributions.H"
#include "include/crossVersionCompatibility.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class porosityCoefficient Declaration
\*---------------------------------------------------------------------------*/

class waveGauges
{
private:

    // Private data

        const fvMesh& mesh_;
        const dictionary& gaugeDict_;


    // Private Member Functions

        void writeVTKFormat
        (
            const word&,
            const pointField&,
            const point&
        );


public:

    //- Runtime type information
    TypeName("waveGauges");


    // Constructors

        //- Construct from components
        waveGauges
        (
            const fvMesh&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        waveGauges(const waveGauges&) = delete;


    //- Destructor
    virtual ~waveGauges()
    {};


    // Member Functions

        // Access
        void evaluate(const word&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const waveGauges&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
