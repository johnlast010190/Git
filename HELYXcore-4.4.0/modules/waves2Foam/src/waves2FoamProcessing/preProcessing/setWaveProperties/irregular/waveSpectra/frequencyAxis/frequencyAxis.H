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
    Foam::frequencyAxis

Description
    An abstract base class for setting the frequencyAxis for irregular wave trains.

SourceFiles
    frequencyAxis.C

Author
    Niels Gjøl Jacobsen, Deltares.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef frequencyAxis_H
#define frequencyAxis_H

#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "preProcessing/setWaveProperties/setWaveProperties.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class frequencyAxis Declaration
\*---------------------------------------------------------------------------*/

class frequencyAxis
{
protected:

    // Protected data

        const Time& rT_;

        dictionary& dict_;

        scalar fp_;
        scalar fl_;
        scalar fu_;


public:

    //- Runtime type information
    TypeName("frequencyAxis");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            frequencyAxis,
            frequencyAxis,
            (
                const Time& rT,
                dictionary& dict
            ),
            (rT, dict)
        );


    // Constructors

        //- Construct from components
        frequencyAxis
        (
            const Time&,
            dictionary&
        );

        //- Disallow default bitwise copy construction
        frequencyAxis(const frequencyAxis&) = delete;


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<frequencyAxis> New
        (
            const Time&,
            dictionary&
        );


    //- Destructor
    virtual ~frequencyAxis();


    // Member functions

        virtual scalarField freqAxis
        (
            const scalarField&,
            const scalarField&,
            const label&
        ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const frequencyAxis&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
