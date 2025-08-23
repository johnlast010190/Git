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
    (c) 2012-2016 OpenFOAM Foundation
    (c) 2023 FOSS GP
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMotionSolvers/sinus/sinusFvMotionSolver.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sinusFvMotionSolver, 0);
    addToRunTimeSelectionTable
    (
        motionSolver,
        sinusFvMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sinusFvMotionSolver::sinusFvMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    componentVelocityMotionSolver(mesh, dict, typeName),
    fvMotionSolver(mesh),
    period_(coeffDict().lookup<scalar>("period")),
    radFrequency_(2*M_PI/period_),
    phase_(coeffDict().lookup<scalar>("phase")),
    amplitude_(coeffDict().lookup<scalar>("amplitude")),
    points0_(fvMesh_.points())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::sinusFvMotionSolver::curPoints() const
{
    const pointField& points = fvMesh_.points();
    tmp<pointField> tcurPoints(new pointField(points));

    scalarField amplitude_x(points0_.component(cmpt_));

    forAll(amplitude_x, i)
    {
        if (amplitude_x[i] < -0.5)
        {
            amplitude_x[i] = (amplitude_x[i]+4.0) / 4.0;
        }
        else if (amplitude_x[i] >= -0.5 && amplitude_x[i] <= 0.5)
        {
            amplitude_x[i] = 1.0;
        }
        else if (amplitude_x[i] > 0.5)
        {
            amplitude_x[i] = 1.0 - ((amplitude_x[i]-0.5)/ 4.0) ;
        }
    }

    tcurPoints.ref().replace
    (
        cmpt_,
        points0_.component(cmpt_)+amplitude_x
       *(amplitude_*Foam::sin(radFrequency_*fvMesh_.time().value() + phase_))
    );

    twoDCorrectPoints(tcurPoints.ref());

    return tcurPoints;
}


// ************************************************************************* //
