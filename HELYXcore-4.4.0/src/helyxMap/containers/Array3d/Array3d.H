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

Description:
    STL-based three-dimensional container for data storage and manipulation.

\*---------------------------------------------------------------------------*/

#ifndef Array3d_H
#define Array3d_H

#include <vector>
#include "containers/Array2d/Array2d.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
struct Tri
{
    Type a;
    Type b;
    Type c;

    Tri(const Type& a0, const Type& b0, const Type& c0)
    {
        a = a0;
        b = b0;
        c = c0;
    }
};


/*---------------------------------------------------------------------------*\
                           Class Array3d Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Array3d
{
public:

    // Public Data

        std::vector<Array2d<Type>> data_;


    // Constructors

        //- Default construct
        Array3d()
        {}


    // Member Functions

        void resize(label n1)
        {
            data_.resize(n1);
        }

        void setSize(label nx, label ny, label nz)
        {
            data_.resize(nx);

            for (label i = 0; i < nx; i++)
            {
                data_[i].resize(ny);
            }

            for (label i = 0; i < nx; i++)
            {
                for (label j = 0; j < ny; j++)
                {
                    data_[i][j].resize(nz);
                }
            }
        }

        void pushBack(Array2d<Type>& data)
        {
            data_.push_back(data);
        }

        label dim1() const
        {
            return data_.size();
        }

        label dim2() const
        {
            return data_[0].size();
        }

        label dim3() const
        {
            return data_[0][0].size();
        }

        void setVal(label i, label j, label k, const Type& val)
        {
            data_[i][j][k] = val;
        }

        Type getVal(label i, label j, label k)
        {
           return data_[i][j][k];
        }

        label size() const
        {
            return data_.size();
        }

        bool empty()
        {
            return data_.empty();
        }


    // Member Operators

        Array2d<Type>& operator[](label i)
        {
            return data_[i];
        }

        const Array2d<Type>& operator[](label i) const
        {
            return data_[i];
        }

        Type& operator()(label i, label j, label k)
        {
            return data_[i][j][k];
        }

        Type operator()(label i, label j, label k) const
        {
            return data_[i][j][k];
        }

        //- Assign a value to the whole array
        inline void operator=(const Type& val);

        void operator+=(const Type& val);

        void operator-=(const Type& val);

        void operator*=(const Type& val);

        void operator/=(const Type& val);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Array3dI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
