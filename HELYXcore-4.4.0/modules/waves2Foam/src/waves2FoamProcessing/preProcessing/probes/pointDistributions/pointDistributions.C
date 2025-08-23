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
    (c) held by original author

\*---------------------------------------------------------------------------*/

#include "pointDistributions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(pointDistributions, 0);
defineRunTimeSelectionTable(pointDistributions, pointDistributions);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


pointDistributions::pointDistributions
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),

    pointDict_( dict )
{
}


autoPtr<pointDistributions> pointDistributions::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    word pd( dict.lookup("pointDistribution") );

    const auto ctor =
        ctorTableLookup
        (
            "point distribution",
            pointDistributionsConstructorTable_(),
            pd
        );
    return autoPtr<pointDistributions>(ctor( mesh, dict));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
