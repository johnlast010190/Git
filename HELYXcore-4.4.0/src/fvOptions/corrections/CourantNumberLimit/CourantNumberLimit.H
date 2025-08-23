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
    (c) 2017-2024 Engys Ltd.
    (c) 2015 OpenCFD Ltd.
    (c) 2015 OpenFOAM Foundation

Class
    Foam::fv::CourantNumberLimit

Group
    grpFvOptionsConstraints

Description
    Constraint to clip velocity extrema in transient calculations

    Velocity is limited based on Courant number. No source or sinks are
    applied, excessive velocity is simply scaled to comply with the
    constraint CFL.

    Constraint described by:

        CourantNumberLimitCoeffs
        {
            cflMax        100;

            // Optional: name of velocity field (default: U)
            //UName       U;

            //Optional flux name
            //phiName phi;

            //Optional phase field name
            //phaseName none;

            //Optional add scaled clipped velocity back into field
            //residualScale 0;

        }


SourceFiles
    CourantNumberLimit.C

\*---------------------------------------------------------------------------*/

#ifndef CourantNumberLimit_H
#define CourantNumberLimit_H

#include "cellSetOption/cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class CourantNumberLimit Declaration
\*---------------------------------------------------------------------------*/

class CourantNumberLimit
:
    public cellSetOption
{

protected:

    // Protected data

        //- Maximum Courant number
        scalar cflMax_;

        //- scale clipped portion of velocity, useful for diagnostics
        scalar residualScale_;

        //- Flux field used to calculate Co number
        word phiName_;

        //- multiphase field for interface CFL limiting
        word phaseName_;

        //- density name
        word rhoName_;


public:

    //- Runtime type information
    TypeName("CourantNumberLimit");


    // Constructors

        //- Construct from components
        CourantNumberLimit
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        CourantNumberLimit(const CourantNumberLimit&) = delete;


    //- Destructor
    virtual ~CourantNumberLimit()
    {}


    // Member Functions

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames);

        //- Clip values
        virtual void correct(volVectorField& U);

        //- Read dictionary
        virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const CourantNumberLimit&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
