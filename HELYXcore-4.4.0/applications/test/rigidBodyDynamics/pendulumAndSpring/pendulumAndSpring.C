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
    (c) 2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Application
    pendulumAndSpring

Description
    Simple 2-DoF pendulum and spring simulation.

\*---------------------------------------------------------------------------*/

#include "rigidBodyMotion/rigidBodyMotion.H"
#include "bodies/masslessBody/masslessBody.H"
#include "bodies/sphere/sphere.H"
#include "joints.H"
#include "restraints/restraint/rigidBodyRestraint.H"
#include "rigidBodyModelState/rigidBodyModelState.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "global/constants/constants.H"

using namespace Foam;
using namespace RBD;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    dictionary pendulumAndSpringDict(IFstream("pendulumAndSpring")());

    // Create the pendulumAndSpring model from dictionary
    rigidBodyMotion pendulumAndSpring(pendulumAndSpringDict);

    label nIter(pendulumAndSpringDict.lookup<label>("nIter"));

    Info<< pendulumAndSpring << endl;
    Info<< "// Joint state " << endl;
    pendulumAndSpring.state().write(Info);

    // Create the joint-space force field
    scalarField tau(pendulumAndSpring.nDoF(), Zero);

    // Create the external body force field
    Field<spatialVector> fx(pendulumAndSpring.nBodies(), Zero);

    OFstream xFile("xVsTime");
    OFstream omegaFile("omegaVsTime");

    // Integrate the motion of the pendulumAndSpring for 100s
    scalar deltaT = 0.01;
    for (scalar t=0; t<20; t+=deltaT)
    {
        pendulumAndSpring.newTime();

        for (label i=0; i<nIter; i++)
        {
            pendulumAndSpring.solve(deltaT, tau, fx);
        }

        // Write the results for graph generation
        xFile << t << " " << pendulumAndSpring.state().q()[0] << endl;
        omegaFile << t << " " << pendulumAndSpring.state().q()[1] << endl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
