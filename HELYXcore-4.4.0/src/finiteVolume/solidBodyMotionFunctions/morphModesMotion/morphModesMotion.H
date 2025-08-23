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
    Foam::solidBodyMotionFunctions::morphModesMotion

Description
    Motion that follows the modes of a morphing mesh. However, in future we
    want to extend this to follow mesh deformation motion as well. To avoid
    user clutter it would be beneficial to move some of the functionality
    to a base class and share it with other deforming motion classes
    (example is the search for a points on the patch).

    Example:
    \verbatim
    motionFunction
    {
        type            morphModesMotion;
        files           ( pointMotionU_rot );
        pointMotionU_rot table
        (
            ( 0 1 )
            ( 0.16667 0 )
            ( 0.33333 1 )
            ( 0.5 0 )
            ( 0.66667 1 )
            ( 0.83333 0 )
            ( 1 1 )
        );
        patches ( motorBike );
        incrementalMotion false;
    }
    \endverbatim

    Note: Motion is still missing proper restart from time0 points, because
    it seems like the modesMotionSolver is missing it as well and this motion
    is intended to work in conjungtion with the modesMotionSolver.

SourceFiles
    morphModesMotion.C

\*---------------------------------------------------------------------------*/

#ifndef morphModesMotion_H
#define morphModesMotion_H

#include "solidBodyMotionFunctions/solidBodyMotionFunction/solidBodyMotionFunction.H"
#include "fields/Fields/primitiveFields.H"
#include "meshes/primitiveShapes/point/point.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "memory/autoPtr/autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                    Class morphModesMotion Declaration
\*---------------------------------------------------------------------------*/

class morphModesMotion
:
    public solidBodyMotionFunction
{
    // Private data

        //- Mesh reference
        const fvMesh& mesh_;

        //- List of names containing the modes
        wordReList fileNames_;

        //- Amplification of each mode
        PtrList<Function1<scalar>> amplifications_;

        //- Point vector values for morphing mesh on picked triangle
        List<vectorField> triangleModes_;

        //- CofR of the morphing mesh
        mutable vector CofR_;

        //- CofR of the morphing mesh
        mutable vector CofR0_;

        //- Triangle points
        mutable vectorField triPoints_;

        //- Triangle points
        mutable vectorField triPoints0_;


public:

    //- Runtime type information
    TypeName("morphModesMotion");


    // Constructors

        //- Construct from components
        morphModesMotion
        (
            const objectRegistry& obr,
            const dictionary& SBMFCoeffs,
            const word& frameName
        );

        //- Disallow default bitwise copy construction
        morphModesMotion(const morphModesMotion&) = delete;

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const override
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new morphModesMotion(time_.db(), SBMFCoeffs_, frameName_)
            );
        }


    //- Destructor
    virtual ~morphModesMotion()
    {}


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation
        (
            const scalar t1,
            const scalar t2
        ) const override;

        //- Return the motion velocity
        virtual vectorTuple velocity() const override
        {
            return vectorTuple(Zero, Zero);
        }

        //- Return the solid-body motion transformation septernion
        labelList findTrianglePoints(const vectorField& points);

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs) override;

        //- Center of rotation
        virtual vector& CofR() const override
        {
            return CofR_;
        }

        //- Center of rotation
        virtual vector& CofR0() const override
        {
            return CofR0_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const morphModesMotion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
