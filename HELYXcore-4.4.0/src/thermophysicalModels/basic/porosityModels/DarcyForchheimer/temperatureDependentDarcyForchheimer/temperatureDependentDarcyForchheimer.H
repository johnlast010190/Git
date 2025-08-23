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
    (c) 2010-2019 Engys Ltd.
    (c) 2024 Engys Ltd.

Description
    temperature dep. Darcy-Forchheimer law porosity model, given by:

        \f[
            S = - (\mu d*d(T)/d + \frac{\rho |U|}{2} f*f(T)/f) U
        \f]

        equivalent to

        \f[
            S = - (\mu d(T) + \frac{\rho |U|}{2} f(T)) U
        \f]

    where
    \vartable
        d(T)/d   | viscous, as Function1
        f(T)/f   | inertial, as Function1
    \endvartable

    If pureUdP_= true, a pure dP vs U relation
    is used

        \f[
            S = - (d(T) + f(T)*U ) U
        \f]

SourceFiles
    temperatureDependentDarcyForchheimer.C
    temperatureDependentDarcyForchheimerTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef temperatureDependentDarcyForchheimer_H
#define temperatureDependentDarcyForchheimer_H

#include "porosityModels/DarcyForchheimer/DarcyForchheimer.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace porosityModels
{

/*---------------------------------------------------------------------------*\
                      Class DarcyForchheimer Declaration
\*---------------------------------------------------------------------------*/

class temperatureDependentDarcyForchheimer
:
    public DarcyForchheimer
{
public:

    //- Runtime type information
    TypeName("temperatureDependentDarcyForchheimer");

    // Constructors

        //- Constructor
        temperatureDependentDarcyForchheimer
        (
            const word& name,
            const word& modelType,
            const objectRegistry& obr,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& cellZoneName
        );

        //- Disallow default bitwise copy construction
        temperatureDependentDarcyForchheimer
        (
            const temperatureDependentDarcyForchheimer&
        ) = delete;

    //- Destructor
    virtual ~temperatureDependentDarcyForchheimer();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const temperatureDependentDarcyForchheimer&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace porosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
