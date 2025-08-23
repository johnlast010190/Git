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
    Foam::interpolationWeights

Description
    Abstract base class for interpolating in 1D

SourceFiles
    interpolationWeights.C
    interpolationWeightsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef interpolationWeights_H
#define interpolationWeights_H

#include "fields/Fields/scalarField/scalarField.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "meshes/primitiveShapes/point/pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class objectRegistry;

/*---------------------------------------------------------------------------*\
                       Class interpolationWeights Declaration
\*---------------------------------------------------------------------------*/

class interpolationWeights
{
protected:

        const scalarField& samples_;

public:

    //- Runtime type information
    TypeName("interpolationWeights");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            interpolationWeights,
            word,
            (
                const scalarField& samples
            ),
            (samples)
        );


    // Constructors

        //- Construct from components
        interpolationWeights(const scalarField& samples);

        //- Disallow default bitwise copy construction
        interpolationWeights(const interpolationWeights&) = delete;


    // Selectors

        //- Return a reference to the selected interpolationWeights
        static autoPtr<interpolationWeights> New
        (
            const word& type,
            const scalarField& samples
        );


    //- Destructor
    virtual ~interpolationWeights();


    // Member Functions

        //- Calculate weights and indices to calculate t from samples.
        //  Returns true if indices changed.
        virtual bool valueWeights
        (
            const scalar t,
            labelList& indices,
            scalarField& weights
        ) const = 0;

        //- Calculate weights and indices to calculate integrand of t1..t2
        //  from samples. Returns true if indices changed.
        virtual bool integrationWeights
        (
            scalar t1,
            scalar t2,
            labelList& indices,
            scalarField& weights
        ) const = 0;

        //- Calculate weights and indices to calculate derivative of t1..t2
        //  from samples. Returns true if indices changed.
        virtual bool derivationWeights
        (
            const scalar x,
            labelList& indices,
            scalarField& weights
        ) const = 0;

        //- Helper: weighted sum
        template<class ListType1, class ListType2>
        static typename outerProduct
        <
            typename ListType1::value_type,
            typename ListType2::value_type
        >::type
        weightedSum(const ListType1& f1, const ListType2& f2);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const interpolationWeights&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolations/interpolationWeights/interpolationWeights/interpolationWeightsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
