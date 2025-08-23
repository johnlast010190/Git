/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) held by original author
    (c) 2024 Engys Ltd.

Class
    Foam::solidBodyMotionFunctions::swakMotion

Description
    Returns a solid body motion according to an evaluation of swak-functions

SourceFiles
    swakMotion.C

\*---------------------------------------------------------------------------*/

#ifndef swakMotion_H
#define swakMotion_H

#include "solidBodyMotionFunctions/solidBodyMotionFunction/solidBodyMotionFunction.H"
#include "fields/Fields/primitiveFields.H"
#include "meshes/primitiveShapes/point/point.H"

#include "CommonValueExpressionDriver.H"
#include "helpers/exprString.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class swakMotion Declaration
\*---------------------------------------------------------------------------*/

class swakMotion
:
    public solidBodyMotionFunction
{
    // Private data

        //- The actual driver for the Expressions
        autoPtr<CommonValueExpressionDriver> driver_;

        bool doTranslation_;

        bool doRotation_;

        bool alphaIsDegrees_;

        exprString translationExpression_;

        exprString axisExpression_;

        exprString originExpression_;

        exprString alphaExpression_;


    // Private Member Functions

        CommonValueExpressionDriver& driver() const;


public:

    //- Runtime type information
    TypeName("swakMotion");


    // Constructors

        //- Construct from components
        swakMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime,
            const word& frameName = "none"
        );

        //- Disallow default bitwise copy construction
        swakMotion(const swakMotion&) = delete;

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const override
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new swakMotion(SBMFCoeffs_, time_, frameName_)
            );
        }


    //- Destructor
    virtual ~swakMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const override;

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation
        (
            const scalar t0,
            const scalar t
        ) const override
        {
            NotImplemented;
        }

        //- Return the solid-body motion velocity
        // transformation septernion (for foam-extend)
        virtual vectorTuple velocity() const override;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs) override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const swakMotion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
