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
    (c) 2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::tabulated6DoFAcceleration

Description
    Tabulated 6DoF acceleration.

    Obtained by interpolating tabulated data for linear acceleration,
    angular velocity and angular acceleration.

SourceFiles
    tabulated6DoFAcceleration.C

\*---------------------------------------------------------------------------*/

#ifndef tabulated6DoFAcceleration_H
#define tabulated6DoFAcceleration_H

#include "fields/Fields/primitiveFields.H"
#include "primitives/Vector2D/Vector2D.H"
#include "db/Time/Time.H"
#include "meshes/primitiveShapes/point/point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class tabulated6DoFAcceleration Declaration
\*---------------------------------------------------------------------------*/

class tabulated6DoFAcceleration
{
    // Private data

        const Time& time_;

        dictionary accelerationCoeffs_;

        //- Time data file name read from dictionary
        fileName timeDataFileName_;

        //- Type used to read in the acceleration "vectors"
        typedef Vector<vector> accelerationVectors;

        //- Field of times
        scalarField times_;

        //- Field of acceleration "vectors"
        Field<accelerationVectors> values_;

        //- optional origin
        point origin_;


public:

    //- Runtime type information
    TypeName("tabulated6DoFAcceleration");


    // Constructors

        //- Construct from components
        tabulated6DoFAcceleration
        (
            const dictionary& accelerationCoeffs,
            const Time& runTime
        );

        //- Disallow default bitwise copy construction
        tabulated6DoFAcceleration(const tabulated6DoFAcceleration&) = delete;


    //- Destructor
    virtual ~tabulated6DoFAcceleration();


    // Member Functions

        //- Return the solid-body accelerations
        virtual Vector<vector> acceleration() const;

        //- Return the origin location
        virtual point origin() const
        {
            return origin_;
        }

        //- Update properties from given dictionary
        virtual bool read(const dictionary& accelerationCoeffs);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const tabulated6DoFAcceleration&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
