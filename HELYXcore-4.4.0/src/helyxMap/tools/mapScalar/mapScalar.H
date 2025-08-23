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
    (c) 2024-2025 Engys Ltd.

Class
    Foam::mapScalar

Description:
    Map a scalar field.

SourceFiles
    mapScalar.C

\*---------------------------------------------------------------------------*/

#ifndef mapScalar_H
#define mapScalar_H

#include "knn/knnInterpolation.H"
#include "referenceFrames/coordinateFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
namespace incompressible
{
    class turbulentConvectiveTemperatureFvPatchScalarField;
}
namespace compressible
{
    class turbulentConvectiveTemperatureFvPatchScalarField;
}

/*---------------------------------------------------------------------------*\
                          Class mapScalar Declaration
\*---------------------------------------------------------------------------*/

class mapScalar
{
    // Private Data

        //- Name of the field data table, defaults to the name of the field.
        word fieldTableName_;

        //- If true, points are in cylindrical coordinates and must be
        //  converted to cartesian
        Switch cylindricalCoords_;

        Switch f4stFormat_;

        autoPtr<knnInterpolation> mapperPtr_;

        //- Need to store sample points for use with frames
        autoPtr<pointField> samplePoints_;

        //- List of boundaryData time directories
        instantList sampleTimes_;

        const scalar timeValue_;

        //- Current starting index in sampleTimes
        label startSampleTime_;

        //- Number of neighbours used in kNN
        label numNbrs_;


public:

    friend class
        incompressible::turbulentConvectiveTemperatureFvPatchScalarField;
    friend class compressible::turbulentConvectiveTemperatureFvPatchScalarField;


    // Constructors

        //- Default construct
        mapScalar();

        //- Construct from the field name
        mapScalar(const word& fieldName);

        mapScalar(const mapScalar& mapsc, const word& fieldName);

        mapScalar
        (
            const word& fieldName,
            const dictionary& dict,
            const scalar timeValue,
            const label neighbours = 5
        );


    // Member Functions

        //- Transform points from local to global frame
        tmp<pointField> makePointsGlobal
        (
            pointField& samplePoints,
            const coordinateFrame* coordFrPtr,
            bool isDefined
        );

        void checkTable
        (
            const Time& time,
            const fvPatch& patch,
            const coordinateFrame* coordFrPtr,
            bool isDefined
        );

        tmp<scalarField> interpolateSampledData
        (
            const Time& time,
            const fvPatch& patch
        );

        void clearInterpolator()
        {
            mapperPtr_.clear();
            startSampleTime_ = -1;
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
