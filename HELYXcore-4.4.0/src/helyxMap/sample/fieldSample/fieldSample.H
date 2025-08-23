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
    (c) 2019-2025 Engys Ltd.

Class
    Foam::fieldSample

Description
    Class designed to sample scalar and vector fields from an unstructured mesh.

SourceFiles
    fieldSample.C

\*---------------------------------------------------------------------------*/

#ifndef fieldSample_H
#define fieldSample_H

#include "sample/gridField/gridField.H"
#include "fields/volFields/volFields.H"
#include "src/knn.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fieldSample Declaration
\*---------------------------------------------------------------------------*/

class fieldSample
:
    public gridField
{
public:

    // Public Data

        KNN bodyMap_;

        KNN internalMap_;

        word bodySurfaceType_;

        word bodySurfaceName_;

        bool byName_;

        bool byNameList_;

        bool byPhysicalType_;

        const volScalarField* field_;

        DynamicList<std::pair<label, label>> bodyFaces_;

        Array3d<scalar> gridField_;


    // Constructors

        //- Default construct
        fieldSample();

        //- Construct from components
        fieldSample
        (
            label nx,
            label ny,
            label nz,
            scalar x0,
            scalar y0,
            scalar z0,
            scalar x1,
            scalar y1,
            scalar z1
        );


    //- Destructor
    virtual ~fieldSample();


    // Member Functions

        void writeSamples(const word& fldName);

        void sampleField
        (
            const volScalarField& fld,
            const fvMesh* mesh,
            bool setLevel = false    // set distance inside solid
        );

        inline void getSolidField
        (
            label i,
            label j,
            label k,
            const std::vector<label>& nb,
            const std::vector<scalar>& dist,
            bool setLevel = false
        );

        //- Get the scalar field at a nodal point by kdTree search
        inline void getGridField
        (
            label i,
            label j,
            label k,
            const std::vector<label>& nb,
            const std::vector<scalar>& dist
        );

        void constructFieldKnn(const fvMesh* mesh);

        void constructBodyKnn(const fvMesh* mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fieldSampleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
