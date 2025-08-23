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
    (c) 2015-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "finiteVolume/ddtSchemes/localEulerDdtScheme/localEulerDdtScheme.H"
#include "finiteVolume/ddtSchemes/boundedDdtScheme/boundedDdtScheme.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::word Foam::fv::localEulerDdt::rDeltaTName("rDeltaT");
Foam::word Foam::fv::localEulerDdt::rDeltaTfName("rDeltaTf");
Foam::word Foam::fv::localEulerDdt::rSubDeltaTName("rSubDeltaT");

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::fv::localEulerDdt::enabled(const fvMesh& mesh)
{
    ITstream& s(mesh.schemes().ddtScheme("default"));
    word firstWord(s);
    return
        firstWord == fv::localEulerDdtScheme<scalar>::typeName
     || (
            firstWord == word(fv::boundedDdtScheme<scalar>::typeName_())
         && word(s) == fv::localEulerDdtScheme<scalar>::typeName
        );
}


const Foam::volScalarField& Foam::fv::localEulerDdt::localRDeltaT
(
    const fvMesh& mesh
)
{
    return mesh.objectRegistry::lookupObject<volScalarField>
    (
        mesh.time().subCycling() ? rSubDeltaTName : rDeltaTName
    );
}


const Foam::surfaceScalarField& Foam::fv::localEulerDdt::localRDeltaTf
(
    const fvMesh& mesh
)
{
    return mesh.objectRegistry::lookupObject<surfaceScalarField>
    (
        rDeltaTfName
    );
}


Foam::tmp<Foam::volScalarField> Foam::fv::localEulerDdt::localRSubDeltaT
(
    const fvMesh& mesh,
    const label nAlphaSubCycles
)
{
    // Should still be registred in object registry for lookups
    return tmp<volScalarField>
    (
        new volScalarField
        (
            rSubDeltaTName,
            nAlphaSubCycles
           *mesh.objectRegistry::lookupObject<volScalarField>(rDeltaTName)
        )
    );
}


// ************************************************************************* //
