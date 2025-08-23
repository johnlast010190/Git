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
    Foam::gridField

Description
    Class designed to create a symmetric scalar field from an unsymmetric field.

SourceFiles
    gridField.C

\*---------------------------------------------------------------------------*/

#ifndef gridField_H
#define gridField_H

#include "containers/Array3d/Array3d.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "meshes/boundBox/boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class gridField Declaration
\*---------------------------------------------------------------------------*/

class gridField
{
protected:

    // Protected Data

        scalar xmin_;
        scalar xmax_;
        scalar ymin_;
        scalar ymax_;
        scalar zmin_;
        scalar zmax_;


public:

    // Static Data Members

        static scalar rhoRef_;
        static vector Uref_;


    // Public Typedefs

        typedef DynamicList<scalar> nodeField;

        typedef DynamicList<vector> nodeVectors;


    // Public Data

        bool containSolidObject_;

        label nx_;
        label ny_;
        label nz_;

        scalar dx_;
        scalar dy_;
        scalar dz_;

        word name_;

        bool interp_;
        bool interp2_;

        wordList scalarSamples_;
        wordList vectorSamples_;

        Array3d<label> inSolids_;

        Array3d<nodeField> gridFields_;

        Array3d<scalar> dists_;

        Array3d<point> gridNodes_;

        Array3d<nodeVectors> gridVectors_;


    // Static Member Functions

        static label inList(const word& key, const wordList& vals);


    // Constructors

        //- Default construct
        gridField();

        gridField
        (
            const word& bbox, // bounding box xmin,ymin,zmin,xmax,ymax,zmax
            label nx,
            label ny,
            label nz = 1,
            bool containObj = true
        );

        gridField
        (
            const boundBox& bbox, // bounding box xmin,ymin,zmin,xmax,ymax,zmax
            label nx,
            label ny,
            label nz = 1,
            bool containObj = true
        );


    //- Destructor
    virtual ~gridField();


    // Member Functions

        inline label scalarId(const word& scname);

        inline label vectorId(const word& vecname);

        void generateNodes();

        inline bool inBlock(scalar x, scalar y, scalar z);

        //- Coordinate of a grid point
        inline point gridPoint(label i, label j, label k);

        void setScalar
        (
            const word& scName,
            label i,
            label j,
            label k,
            scalar val
        )
        {
            label isc = scalarId(scName);

            gridFields_(i, j, k)[isc] = val;
        }

        void setVector
        (
            const word& vecName,
            label i,
            label j,
            label k,
            point val
        )
        {
            label iv = vectorId(vecName);

            gridVectors_(i, j, k)[iv] = val;
        }

        inline scalar p(label i, label j, label k, label ip);

        inline point U(label i, label j, label k, label iv);

        inline scalar getScalar(word scName, scalar x, scalar y, scalar z);

        inline point getVector(word vecName, scalar x, scalar y, scalar z);

        scalar xmin()
        {
            return xmin_;
        }

        scalar xmax()
        {
            return xmax_;
        }

        word name()
        {
            return name_;
        }

        label nzpt()
        {
            if (nz_ == 1)
            {
                return nz_;
            }
            else
            {
                return nz_ + 1;
            }
        }

        scalar ymin()
        {
            return ymin_;
        }

        scalar ymax()
        {
            return ymax_;
        }

        scalar zmin()
        {
            return zmin_;
        }

        scalar zmax()
        {
            return zmax_;
        }

        label nNodes()
        {
            if (nz_ > 1)
            {
                return (nx_ + 1)*(ny_ + 1)*(nz_ + 1);
            }
            else
            {
                return (nx_ + 1)*(ny_ + 1);
            }
        }

        label dimension()
        {
            if (nz_ == 1)
            {
                return 2;
            }
            else
            {
                return 3;
            }
        }

        scalar xbar(scalar x)
        {
            return (x - xmin())/(xmax() - xmin());
        }

        scalar ybar(scalar y)
        {
            return (y - ymin())/(ymax() - ymin());
        }

        scalar zbar(scalar z)
        {
            return (z - zmin())/(zmax() - zmin());
        }

        label cell_ijk(label i, label j, label k)
        {
            if (nz_ == 1)
            {
                return i + j*nx_;
            }
            else
            {
                return i + j*nx_ + k*nx_*ny_;
            }
        }

        label ijk(label i, label j, label k)
        {
            if (nz_ == 1)
            {
                return i + j*(nx_ + 1);
            }
            else
            {
                return i + j*(nx_ + 1) + k*(nx_ + 1)*(ny_ + 1);
            }
        }

        void getBbox(const word& val);

        void getBbox(const boundBox& val);

        void setDelta();

        void setSampleNames
        (
            const HashSet<word>& scalars,
            const HashSet<word>& vectors
        );

        void setSampleNames(const word& scalars, const word& vectors = "none");

        void setDivision(label nx, label ny, label nz = 1);

        inline void getGridId
        (
            label& i,
            label& j,
            label& k,
            scalar x,
            scalar y,
            scalar z
        );

        label size()
        {
            return (nx_ + 1)*(ny_ + 1)*nzpt();
        }

        void setFieldSize();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "gridFieldI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
