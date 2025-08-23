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
    Foam::rawAlphaProbes

Description
    Reads, concatenates and output raw alpha1 probe data.

SourceFiles
    rawAlphaProbes.C

Author
    Bjarne Jensen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef rawAlphaProbes_H
#define rawAlphaProbes_H

#include "postProcessing/postProcessingWaves/postProcessingWaves.H"
#include "include/crossVersionCompatibility.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class rawAlphaProbes
:
    public postProcessingWaves
{
protected:

    // Protected data

        word inputDir_;

        wordList timeDirs_;

        Switch removeDuplicate_;


private:

    // Private Member Functions

        void resizeFields
        (
            List<std::pair<scalar, label>>&,
            List<scalarField>&,
            label
        );

        void writeRawData
        (
            const List<std::pair<scalar, label>>&,
            const scalarField&,
            const scalarField&,
            const scalarField&,
            const List<scalarField>&
        );


public:

    //- Runtime type information
    TypeName("rawAlphaProbes");


    // Constructors

        //- Construct from components
        rawAlphaProbes
        (
            const Time&,
            const dictionary&,
            const word&
        );

        //- Disallow default bitwise copy construction
        rawAlphaProbes(const rawAlphaProbes&) = delete;


    //- Destructor
    virtual ~rawAlphaProbes();


    // Member Functions

        virtual void evaluate();

        void readAlphaProbesData
        (
            List<std::pair<scalar, label>>& timeLabel,
            scalarField& x,
            scalarField& y,
            scalarField& z,
            List<scalarField>& alphas
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const rawAlphaProbes&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
