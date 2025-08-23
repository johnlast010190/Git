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
    (c) 2017-2024 Engys Ltd.

Namespace
    Foam::geometryTransformations

Description
    Namespace for solid-body motions

Class
    Foam::geometryTransformation

Description
    Base class for user defined coordinate transformation - deformation

SourceFiles
    geometryTransformation.C
    dynamicFvMeshNew.C

\*---------------------------------------------------------------------------*/

#ifndef geometryTransformation_H
#define geometryTransformation_H

#include "meshes/primitiveShapes/point/pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class geometryTransformation Declaration
\*---------------------------------------------------------------------------*/

class geometryTransformation
{
protected:

    // Protected data

        //- user input from dynamicGIBFvMesh
        const dictionary& geoTransCoeffs_;


public:

    //- Runtime type information
    TypeName("geometryTransformation");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            geometryTransformation,
            dictionary,
            (const dictionary& SBMFCoeffs),
            (SBMFCoeffs)
        );


    // Constructors

        //- Construct from the SBMFCoeffs dictionary and Time
        geometryTransformation(const dictionary&);

        //- Disallow default bitwise copy construction
        geometryTransformation(const geometryTransformation&) = delete;


    // Selectors

        //- Select constructed from the SBMFCoeffs dictionary and Time
        static autoPtr<geometryTransformation> New(const dictionary&);


    //- Destructor
    virtual ~geometryTransformation(){};


    // Member Functions

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs) = 0;

        //- transformation of the pointField
        virtual tmp<pointField> transformPoints
        (
            const pointField&
        ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const geometryTransformation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
