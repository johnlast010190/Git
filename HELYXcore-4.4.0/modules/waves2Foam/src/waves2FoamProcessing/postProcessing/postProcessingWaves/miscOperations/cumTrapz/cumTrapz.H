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
    Foam::cumTrapz

Description
    Perform a cumulative trapezoidal integration of the data set.

SourceFiles
    cumTrapz.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef cumTrapz_H
#define cumTrapz_H

#include "postProcessing/postProcessingWaves/postProcessingWaves.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class cumTrapz
:
    public postProcessingWaves
{
protected:

    // Protected data

        IOdictionary dataDict_;

        labelList indices_;


private:

    // Private Member Functions

        void evaluateScalar();

        void evaluateVector();

        void writeScalar
        (
            const List<scalarField>&
        );

        void writeVector
        (
            const List<vectorField>&
        );


public:

    //- Runtime type information
    TypeName("cumTrapz");


    // Constructors

        //- Construct from components
        cumTrapz
        (
            const Time&,
            const dictionary&,
            const word&
        );

        //- Disallow default bitwise copy construction
        cumTrapz(const cumTrapz&) = delete;


    //- Destructor
    virtual ~cumTrapz();


    // Member Functions

        virtual void evaluate();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cumTrapz&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
