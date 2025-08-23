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
    Foam::postProcessingWaves

Description
    Abstract base class for all sorts of different post-processing of
    logged/probed data, e.g. ensemble average, spectral analysis, reflection
    analysis, etc.

SourceFiles
    postProcessingWaves.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef postProcessingWaves_H
#define postProcessingWaves_H

#include "fvMesh/fvMesh.H"
#include "cfdTools/general/include/fvCFD.H"
#include "fields/volFields/volFields.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "primitives/complex/complex.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class porosityCoefficient Declaration
\*---------------------------------------------------------------------------*/

class postProcessingWaves
{

protected:

    // Protected data

        const Time& rT_;
        const dictionary& actionProperties_;

        word callName_;
        word actionType_;
        fileName addDir_;
        fileName directDir_;


    // Protected member functions

        void getTimeDirs
        (
            const word&,
            wordList&
        );

        void writeNameDict
        (
            const scalar& dt,
            const wordList& names
        );

        void writeXYZDict
        (
            const scalar&,
            const scalarField&,
            const scalarField&,
            const scalarField&
        );

        void writeIOScalarField
        (
            const scalarField& field,
            const word& name
        );

        void writeIOVectorField
        (
            const vectorField& field,
            const word& name
        );

        scalarField readIOScalarField
        (
            const word& name
        );

        vectorField readIOVectorField
        (
            const word& name
        );

        word dataType();

        scalar readDeltaT
        (
            const dictionary&
        );

        void readIndices
        (
            const dictionary&,
            labelList&
        );

        List<scalarField> readScalarFields
        (
            const labelList& indices
        );

        List<vectorField> readVectorFields
        (
            const labelList& indices
        );

        scalarField equidistantTime
        (
            const List<std::pair<scalar, label >>& timeLabel,
            const dictionary&
        );

        void interpolationWeights
        (
            const List<std::pair<scalar, label>>&,
            const scalarField&,
            scalarField&,
            labelList&,
            labelList&
        );


public:

    //- Runtime type information
    TypeName("postProcessingWaves");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            postProcessingWaves,
            postProcessingWaves,
            (
                const Time& rT,
                const dictionary& actionProp,
                const word& action
            ),
            ( rT, actionProp, action )
        );


    // Constructors

        //- Construct from components
        postProcessingWaves
        (
            const Time&,
            const dictionary&,
            const word&
        );

        //- Disallow default bitwise copy construction
        postProcessingWaves(const postProcessingWaves&) = delete;


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<postProcessingWaves> New
        (
            const Time&,
            const dictionary&,
            const word&
        );


    //- Destructor
    virtual ~postProcessingWaves();


    // Member Functions

        // Access
        virtual void evaluate() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const postProcessingWaves&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// The GCC compiler do not allow these to be overloaded for some reason,
// hence A and B
bool pairSortA( std::pair<scalar, label>, std::pair<scalar, label>);

bool paitSortB( std::pair<scalar, word>, std::pair<scalar, word>);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
