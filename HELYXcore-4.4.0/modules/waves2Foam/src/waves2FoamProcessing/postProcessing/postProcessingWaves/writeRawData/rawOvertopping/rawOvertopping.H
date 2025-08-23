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
    Foam::rawOvertopping

Description
    Reads, concatenates and output raw overtopping data

SourceFiles
    rawOvertopping.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef rawOvertopping_H
#define rawOvertopping_H

#include "postProcessing/postProcessingWaves/postProcessingWaves.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class rawOvertopping
:
    public postProcessingWaves
{
protected:

    // Protected data

        word inputDir_;

        wordList timeDirs_;

        scalar scaleFlux_;

        Switch removeDuplicate_;


private:

    // Private Member Functions

        void resizeFields
        (
            List<std::pair<scalar, label>>&,
            List<vectorField>&,
            label
        );

        void writeRawData
        (
            const List<std::pair<scalar, label>>&,
            const wordList&,
            const List<vectorField>&
        );


public:

    //- Runtime type information
    TypeName("rawOvertopping");


    // Constructors

        //- Construct from components
        rawOvertopping
        (
            const Time&,
            const dictionary&,
            const word&
        );

        //- Disallow default bitwise copy construction
        rawOvertopping(const rawOvertopping&) = delete;


    //- Destructor
    virtual ~rawOvertopping();


    // Member Functions

        virtual void evaluate();

        void readOvertoppingData
        (
            List<std::pair<scalar, label>>&,
            wordList&,
            List<vectorField>&
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const rawOvertopping&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
