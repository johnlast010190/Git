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
    Foam::rawVelocityProbes

Description
    Reads, concatenates and output the raw velocity signals.

SourceFiles
    rawVelocityProbes.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef rawVelocityProbes_H
#define rawVelocityProbes_H

#include "postProcessing/postProcessingWaves/postProcessingWaves.H"
#include "coordinate/rotation/coordinateRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class rawVelocityProbes
:
    public postProcessingWaves
{
protected:

    // Protected data

        word inputDir_;

        wordList timeDirs_;

        Switch removeDuplicate_;

        tensor R_;


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
            const scalarField&,
            const scalarField&,
            const scalarField&,
            const List<vectorField>&
        );


public:

    //- Runtime type information
    TypeName("rawVelocityProbes");


    // Constructors

        //- Construct from components
        rawVelocityProbes
        (
            const Time&,
            const dictionary&,
            const word&
        );

        //- Disallow default bitwise copy construction
        rawVelocityProbes(const rawVelocityProbes&) = delete;


    //- Destructor
    virtual ~rawVelocityProbes();


    // Member Functions

        virtual void evaluate();

        void readVelocityProbeData
        (
            List<std::pair<scalar, label>>& timeLabel,
            scalarField& x,
            scalarField& y,
            scalarField& z,
            List<vectorField>& etas
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const rawVelocityProbes&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
