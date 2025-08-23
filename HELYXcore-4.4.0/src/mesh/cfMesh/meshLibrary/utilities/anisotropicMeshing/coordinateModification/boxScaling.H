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
    (c) Creative Fields, Ltd.
    (c) 2024 Engys Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Class
    boxScaling

Description
    Checks if a box is contained inside the box object

SourceFiles
    boxScaling.C

\*---------------------------------------------------------------------------*/

#ifndef boxScaling_H
#define boxScaling_H

#include "utilities/anisotropicMeshing/coordinateModification/coordinateModification.H"
#include "meshes/primitiveShapes/point/point.H"
#include "db/typeInfo/typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class boxScaling Declaration
\*---------------------------------------------------------------------------*/

class boxScaling
:
    public coordinateModification
{
    // Private data

        //- Centre of the box
        point centre_;

        //- Length of box sides
        vector lengthVec_;

        //- Scaling factors in all directions
        vector scaleVec_;

        //- Min point
        point pMin_;

        //- Max point
        point pMax_;


    // Private member functions

        //- calculate bounding box points
        void calculateBndBox();


public:

    //- Runtime type information
    TypeName("box");


    // Constructors

        //- Null construct
        boxScaling();

        //- Construct from name, cell size, centre and sizes and scaling factors
        boxScaling
        (
            const word& name,
            const point& centre,
            const scalar lengthX,
            const scalar lengthY,
            const scalar lengthZ,
            const scalar scaleX = 1.0,
            const scalar scaleY = 1.0,
            const scalar scaleZ = 1.0
        );

        //- Construct from dictionary
        boxScaling(const word& name, const dictionary& dict);

        //- Construct and return a clone
        virtual autoPtr<coordinateModification> clone
        (
            const boxScaling& bs
        ) const
        {
            return autoPtr<coordinateModification>
            (
                new boxScaling
                (
                    bs.name(),
                    bs.centre_,
                    bs.lengthVec_.x(),
                    bs.lengthVec_.y(),
                    bs.lengthVec_.z(),
                    bs.scaleVec_.x(),
                    bs.scaleVec_.y(),
                    bs.scaleVec_.z()
                )
            );
        }


    // Member Functions

        //- Return the centre of the box
        virtual point origin() const;

        //- Translate the object to the modified coordinates
        //  this is needed for backward transformation
        virtual void translateAndModifyObject(const vector&);

        //- Calculate the displacement vector for box scaling
        virtual vector displacement(const point&) const;

        //- Calculate the displacement vector for box scaling
        virtual vector backwardDisplacement(const point&) const;

        //- Can this modification object be combined with other ones
        virtual bool combiningPossible() const;

        //- Return that "bounding planes" of the scaling region for
        //  the given object
        virtual void boundingPlanes(PtrList<plane>&) const;

    //- Return as dictionary of entries

        dictionary dict(bool ignoreType = false) const;


      // Write

        //- Write
        void write(Ostream&) const;

        //- Write dictionary
        void writeDict(Ostream&, bool subDict = true) const;


    // Member Operators

        //- Assign from dictionary
        void operator=(const dictionary&);


    // IOstream Operators

        Ostream& operator<<(Ostream&) const;
        friend Ostream& operator<<(Ostream&, const boxScaling&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
