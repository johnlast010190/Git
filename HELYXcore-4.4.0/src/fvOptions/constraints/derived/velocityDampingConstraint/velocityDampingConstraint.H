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
    (c) 2015 OpenCFD Ltd.
    (c) 2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::fv::velocityDampingConstraint

Group
    grpFvOptionsConstraints

Description
    Constraint for velocity to dampen velocity fluctuations in
    steady simulations

    In case of velocity exceeding limit applies
    a source sink to remove the excess velocity by adding a term
    that drives the velocity to 0 (note that we cannot use the old
    term velocity since it most likely is already excessive). This
    additional constraint is scaled with (U-UMax)/dt
    where dt is a local time scale (= velocity/cellsize) so it switches off
    if the velocity is below UMax.

    Constraint described by:

        velocityDampingConstraintCoeffs
        {
            UMax        100;

            // Optional: name of velocity field (default: U)
            //UName       U;
        }


SourceFiles
    velocityDampingConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef velocityDampingConstraint_H
#define velocityDampingConstraint_H

#include "cellSetOption/cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class velocityDampingConstraint Declaration
\*---------------------------------------------------------------------------*/

class velocityDampingConstraint
:
    public cellSetOption
{
protected:

    // Protected data

        //- Maximum velocity
        scalar UMax_;


    // Protected Member Functions

        void addDamping(fvMatrix<vector>& eqn);


public:

    //- Runtime type information
    TypeName("velocityDampingConstraint");


    // Constructors

        //- Construct from components
        velocityDampingConstraint
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        velocityDampingConstraint(const velocityDampingConstraint&) = delete;


    //- Destructor
    virtual ~velocityDampingConstraint()
    {}


    // Member Functions

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames);


        // Set values directly

            //- Constrain vector matrix
            virtual void constrain(fvMatrix<vector>& eqn, const label fieldI);


        // I-O

            //- Write data
            virtual void writeData(Ostream&) const;

            //- Read dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const velocityDampingConstraint&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
