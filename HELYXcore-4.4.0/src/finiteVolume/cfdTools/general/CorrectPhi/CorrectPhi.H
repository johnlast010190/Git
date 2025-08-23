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
    (c) 2015-2016 OpenFOAM Foundation

Class
    Foam::CorrectPhi

Description
    Flux correction functions to ensure continuity.

    Required during start-up, restart, mesh-motion etc. when non-conservative
    fluxes may adversely affect the prediction-part of the solution algorithm
    (the part before the first pressure solution which would ensure continuity).
    This is particularly important for VoF and other multi-phase solver in
    which non-conservative fluxes cause unboundedness of the phase-fraction.

SourceFiles
    CorrectPhi.C

\*---------------------------------------------------------------------------*/

#ifndef CorrectPhi_H
#define CorrectPhi_H

#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    class pressureControl;
    class solutionControl;

    //- If the mesh is moving correct the velocity BCs on the moving walls to
    //  ensure the corrected fluxes and velocity are consistent
    void correctUphiBCs
    (
        volVectorField& U,
        surfaceScalarField& phi
    );

    //- If the mesh is moving correct the velocity BCs on the moving walls to
    //  ensure the corrected fluxes and velocity are consistent
    void correctUphiBCs
    (
        const volScalarField& rho,
        volVectorField& U,
        surfaceScalarField& phi
    );

    //- For legacy solvers only! To be removed!
    template<class RAUfType, class DivUType>
    void CorrectPhi
    (
        volVectorField& U,
        surfaceScalarField& phi,
        const volScalarField& p,
        const RAUfType& rAUf,
        const DivUType& divU,
        solutionControl& ximple,
        const label& refCell=0,
        const word& solver=word::null
    );

    template<class RAUfType, class DivUType>
    void CorrectPhi
    (
        volVectorField& U,
        surfaceScalarField& phi,
        const volScalarField& p,
        const RAUfType& rAUf,
        const DivUType& divU,
        solutionControl& ximple,
        const pressureControl& pCtrl,
        const word& solver=word::null
    );

    template<class RAUfType, class DivRhoUType>
    void CorrectPhi
    (
        volVectorField& U,
        surfaceScalarField& phi,
        const volScalarField& p,
        const volScalarField& rho,
        const volScalarField& psi,
        const RAUfType& rAUf,
        const DivRhoUType& divRhoU,
        solutionControl& ximple,
        const word& solver=word::null,
        const bool& needReference=false,
        const label& refCell=0
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cfdTools/general/CorrectPhi/CorrectPhi.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
