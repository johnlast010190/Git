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
    reflectionAnalysis2DLS

Description
    Performs a reflection analysis for a regular signal. It is based on a
    least-squares approach and the wave number is given by the user.

SourceFiles
    reflectionAnalysis2DLS.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef reflectionAnalysis2DLS_H
#define reflectionAnalysis2DLS_H

#include "postProcessing/postProcessingWaves/postProcessingWaves.H"
#include "postProcessing/postProcessingWaves/spectralAnalysis/spectralMethods/spectralMethodsLeastSquaresBased.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class reflectionAnalysis2DLS
:
    public postProcessingWaves
{
protected:

    // Protected data

        IOdictionary dataDict_;

        labelList indices_;

        label N_;

        scalar period_;

        vector waveNumber_;

        pointField x_;


private:

    // Private Member Functions

        void write
        (
            const scalarField&,
            const scalarField&
        );


public:

    //- Runtime type information
    TypeName("reflectionAnalysis2DLS");


    // Constructors

        //- Construct from components
        reflectionAnalysis2DLS
        (
            const Time&,
            const dictionary&,
            const word&
        );

        //- Disallow default bitwise copy construction
        reflectionAnalysis2DLS(const reflectionAnalysis2DLS&) = delete;


    //- Destructor
    virtual ~reflectionAnalysis2DLS();


    // Member Functions

        virtual void evaluate();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const reflectionAnalysis2DLS&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
