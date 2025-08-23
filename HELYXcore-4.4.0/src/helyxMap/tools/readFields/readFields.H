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
    Foam::readFields

Description:
    Base class for mapping two flow fields with different meshes.

SourceFiles
    readFields.C

\*---------------------------------------------------------------------------*/

#ifndef readFields_H
#define readFields_H

#include <map>
#include "src/knn.h"
#include "meshes/smartBoundBox/smartBoundBox.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Scalar overload of the 'farray2d' type
typedef std::vector<std::vector<scalar>> farray2d;

// Forward declaration of classes
class helyxMap;

/*---------------------------------------------------------------------------*\
                          Class readFields Declaration
\*---------------------------------------------------------------------------*/

class readFields
{
    // Private Member Functions

        //- Whether to exclude a given patch type from a boundary mesh operation
        bool excludePatch(const label patchi) const;


protected:

    // Protected Data

        autoPtr<smartBoundBox> bBox_;

        autoPtr<boundBox> gBox_;

        DynamicList<point> xyz_;

        word type_;


public:

    // Public Typedefs

        typedef unsigned short ushort;

        typedef std::map<word, DynamicList<scalar>> ScalarFields;
        typedef std::map<word, DynamicList<point>> VectorFields;


    // Public Data

        ScalarFields scalarFields_;

        VectorFields vectorFields_;

        word mapTime_;

        const fvMesh* mesh_;

        const Time* runTime_;

        helyxMap* fieldMap_;

        bool norm_;

        std::vector<KNN> kdTrees_;

        KNN boundaryMap_;

        KNN internalMap_;

        scalar alphaMax_;

        //- Map from kdTree ID to cell ID
        std::map<label, label> cellMap_;

        label nCells_;

        Array2d<label> alphaMap_;

        DynamicList<scalar> wDists_;

        scalar deltaY_;

        scalar xmin_;
        scalar xmax_;
        scalar ymin_;
        scalar ymax_;
        scalar zmin_;
        scalar zmax_;

        HashSet<word> scFields_;

        HashSet<word> vecFields_;

        scalar rhoRef_;

        scalar UrotDegreeFromSource_;

        vector Uref_;

        bool interp_;

        scalar maxWallDist_;


    // Constructors

        //- Default construct
        readFields();

        //- Construct from components
        readFields
        (
            const fvMesh* mesh,
            const Time* runTime,
            helyxMap* map,
            const word& type = "source"
        );


    //- Destructor
    virtual ~readFields();


    // Member Functions

        DynamicList<point>& xyz()
        {
            return xyz_;
        }

        const DynamicList<point>& xyz() const
        {
            return xyz_;
        }

        word myType()
        {
            return type_;
        }

        void setInputs();

        void setNorm(bool norm)
        {
            norm_ = norm;
        }

        void constructBoundBox(scalar, const word, label);

        void transformPoints(const boundBox& tbox, DynamicList<point>& xyz);

        const smartBoundBox& bBox() const
        {
            return bBox_();
        }

        const boundBox& gBox() const
        {
            return gBox_();
        }

        point pbar(const point& pt)
        {
            point pp(xbar(pt.x()), ybar(pt.y()), zbar(pt.z()));

            return pp;
        }

        scalar U0()
        {
            return mag(Uref_);
        }

        scalar pRef()
        {
            return U0()*U0()*rhoRef_;
        }

        scalar scale(const word& fldName);

        scalar xbar(scalar x)
        {
            if (norm_)
            {
                return (x - xmin_)/(xmax_ - xmin_);
            }
            else
            {
                return x;
            }
        }

        scalar ybar(scalar y)
        {
            if (norm_)
            {
                return (y - ymin_)/(ymax_ - ymin_);
            }
            else
            {
                return y;
            }
        }

        scalar zbar(scalar z)
        {
            if (norm_)
            {
                return (z - zmin_)/(zmax_ - zmin_);
            }
            else
            {
                return z;
            }
        }

        label bfaceIndex(label j)
        {
            return nCells_ + j;
        }

        inline label cellIndex(scalar alphai, label kdTreeId);

        void constructBoundaryKnn();

        void constructInternalKnn(scalar cutoff = 0.8);

        label numBoundFaces();

        void getBoundBox(const boundBox& box);

        void createFields(const word& timeName = "");

        void buildKdTrees
        (
            const volScalarField& y,
            scalar cutoff = 0.8,
            const word& excludeBndType = "none"
        );

        label inGroup(label i)
        {
            return label(wDists_[i]/deltaY_);
        }

        void setWallDist(const volScalarField& y, scalar maxWDist);

        volScalarField& getScalarField(const word& name)
        {
            return mesh_->lookupObjectRef<volScalarField>(name);
        }

        surfaceScalarField& getFlux(const word& name)
        {
            return mesh_->lookupObjectRef<surfaceScalarField>(name);
        }

        const surfaceScalarField& getFlux(const word& name) const
        {
            return mesh_->lookupObjectRef<surfaceScalarField>(name);
        }

        scalar getScalarField(const word& name, label i)
        {
            return mesh_->lookupObjectRef<volScalarField>(name)
                .primitiveField()[i];
        }

        const volScalarField& getScalarField(const word& name) const
        {
            return mesh_->lookupObjectRef<volScalarField>(name);
        }

        volVectorField& getVectorField(const word& name)
        {
            return mesh_->lookupObjectRef<volVectorField>(name);
        }

        const volVectorField& getVectorField(const word& name) const
        {
            return mesh_->lookupObjectRef<volVectorField>(name);
        }

        const point& getVectorField(const word& name, label i) const
        {
            return mesh_->lookupObjectRef<volVectorField>(name)
                .primitiveField()[i];
        }

        //- Store source fields
        void storeFields();

        void storeFields(const boundBox& sampleBox);

        //- Search is always based on relative coordinates
        inline void searchNearest
        (
            std::vector<label>& nb,
            std::vector<scalar>& dist,
            scalar alphai,
            scalar x,
            scalar y,
            scalar z
        );

        inline void searchNearestBFace
        (
            std::vector<label>& nb,
            std::vector<scalar>& dist,
            scalar x,
            scalar y,
            scalar z
        );

        inline bool inDomain(const point& pt);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "readFieldsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
