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
    Foam::Array2d

Description:
    STL-based two-dimensional container for data storage and manipulation.

SourceFiles
    Array2d.C

\*---------------------------------------------------------------------------*/

#ifndef Array2d_H
#define Array2d_H

#include <vector>
#include "primitives/ints/label/label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Array2d Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Array2d
{
public:

    // Public Data

        std::vector<std::vector<Type>> data_;


    // Constructors

        //- Default construct
        Array2d()
        {}

        //- Construct from components
        Array2d(label m, label n)
        {
            resize(m, n);
        }


    // Member Functions

        void setSize(label nx, label ny)
        {
            data_.resize(nx);

            for (label i = 0; i < nx; i++)
            {
                data_[i].resize(ny);
            }
        }

        void resize(label m)
        {
            data_.resize(m);
        }

        unsigned size() const
        {
            return data_.size();
        }

        label dim1() const
        {
            return data_.size();
        }

        label dim2() const
        {
            return data_[0].size();
        }

        void resize(label m, label n)
        {
            data_.resize(m);

            for (label i = 0; i < m; i++)
            {
                data_[i].resize(n);
            }
        }

        void setval(label i, label j, const Type& val)
        {
            data_[i][j] = val;
        }

        bool empty() const
        {
            return data_.empty();
        }


    // Member Operators

        Type& operator()(label i, label j)
        {
            return data_[i][j];
        }

        std::vector<Type>& operator[](label i)
        {
            return data_[i];
        }

        const std::vector<Type>&operator[](label i) const
        {
            return data_[i];
        }

        Type operator()(label i, label j) const
        {
            return data_[i][j];
        }

        void operator=(Type val)
        {
            for (unsigned i = 0; i < data_.size(); i++)
            {
                for (unsigned j = 0; j < data_[i].size(); j++)
                {
                    data_[i][j] = val;
                }
            }
        }

        void operator+=(const Type& val);

        void operator-=(const Type& val);

        void operator*=(const Type& val);

        void operator/=(const Type& val);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
