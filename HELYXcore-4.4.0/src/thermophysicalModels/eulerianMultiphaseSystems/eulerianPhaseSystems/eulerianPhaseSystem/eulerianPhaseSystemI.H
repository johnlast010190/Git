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
    (c) 2014-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::fvMesh& Foam::eulerianPhaseSystem::mesh() const
{
    return mesh_;
}


inline const Foam::eulerianPhaseSystem::phaseModelList&
Foam::eulerianPhaseSystem::phases() const
{
    return phaseModels_;
}


inline Foam::eulerianPhaseSystem::phaseModelList&
Foam::eulerianPhaseSystem::phases()
{
    return phaseModels_;
}


inline const Foam::eulerianPhaseSystem::phasePairTable&
Foam::eulerianPhaseSystem::phasePairs() const
{
    return phasePairs_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::eulerianPhaseSystem::alphas() const
{
    return thermo_.alphas();
}


inline Foam::PtrList<Foam::volScalarField>&
Foam::eulerianPhaseSystem::alphas()
{
    return thermo_.alphas();
}


inline const Foam::surfaceScalarField& Foam::eulerianPhaseSystem::phiv() const
{
    return phiv_;
}


inline Foam::surfaceScalarField& Foam::eulerianPhaseSystem::phiv()
{
    return phiv_;
}


inline const Foam::volScalarField& Foam::eulerianPhaseSystem::dpdt() const
{
    return dpdt_;
}


inline Foam::volScalarField& Foam::eulerianPhaseSystem::dpdt()
{
    return dpdt_;
}


inline Foam::fv::options& Foam::eulerianPhaseSystem::fvOptions() const
{
    return fv::options::New(mesh_);
}


// ************************************************************************* //
