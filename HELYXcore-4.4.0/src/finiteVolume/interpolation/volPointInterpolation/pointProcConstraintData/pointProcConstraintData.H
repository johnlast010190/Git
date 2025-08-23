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
    (c) 2024 Engys Ltd.

Class

Description

SourceFiles
    pointProcConstraintData.C

\*---------------------------------------------------------------------------*/

#ifndef pointProcConstraintData_H
#define pointProcConstraintData_H

#include "db/typeInfo/typeInfo.H"
#include "meshes/primitiveShapes/point/point.H"
#include "primitives/transform/transform.H"
#include "primitives/transform/transformer/transformer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
template<class Type>
class pointProcConstraintData;

template<class Type>
Istream& operator>>(Istream&, pointProcConstraintData<Type>&);

template<class Type>
Ostream& operator<<(Ostream&, const pointProcConstraintData<Type>&);

template<class Type>
pointProcConstraintData<Type> max
(
    const pointProcConstraintData<Type>&,
    const pointProcConstraintData<Type>&
);

/*---------------------------------------------------------------------------*\
                          Class pointProcConstraintData Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class pointProcConstraintData
{
    // Private data

    // Is proc point?
    bool procPoint_;

    // Is proc point?
    bool procBoundaryPoint_;

    // Does the point fix value from other BCs?
    bool fixesValue_;

    // Point value if it is fixed
    Type fixedValue_;

public:


    // Constructors

        //- Construct from components
        pointProcConstraintData
        (
            const bool&,
            const bool&,
            const bool&,
            const Type&
        );
        pointProcConstraintData();

    // Member Functions


        //- update

        void checkParallelBoundaryIntersections();

        //- Access

        inline bool& procPoint()
        {
            return procPoint_;
        }

        inline const bool& procPoint() const
        {
            return procPoint_;
        }

        inline bool& procBoundaryPoint()
        {
            return procBoundaryPoint_;
        }

        inline const bool& procBoundaryPoint() const
        {
            return procBoundaryPoint_;
        }


        inline bool& fixesValue()
        {
            return fixesValue_;
        }

        inline const bool& fixesValue() const
        {
            return fixesValue_;
        }

        inline Type& fixedValue()
        {
            return fixedValue_;
        }

        inline const Type& fixedValue() const
        {
            return fixedValue_;
        }


    // Member Operators

         bool operator==(const pointProcConstraintData&) const;
         bool operator!=(const pointProcConstraintData&) const;
         bool operator>(const pointProcConstraintData&) const;
         bool operator<(const pointProcConstraintData&) const;


    // IOstream Operators

        friend Istream& operator>> <Type>(Istream&, pointProcConstraintData&);
        friend Ostream& operator<< <Type>(Ostream&, const pointProcConstraintData&);

};

#define DefineContiguousPointProcConstraintDataType(Type, nullArg)             \
    template<>                                                                 \
    constexpr bool contiguous<pointProcConstraintData<Type>>()                 \
    {                                                                          \
        return true;                                                           \
    }

FOR_ALL_FIELD_TYPES(DefineContiguousPointProcConstraintDataType);

#undef DefineContiguousPointProcConstraintDataType


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
  //
#ifdef NoRepository
#include "pointProcConstraintDataI.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
