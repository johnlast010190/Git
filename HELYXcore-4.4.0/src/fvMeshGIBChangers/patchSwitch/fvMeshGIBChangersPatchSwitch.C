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
    (c) 2024-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "patchSwitch/fvMeshGIBChangersPatchSwitch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshGIBChangers
{
    defineTypeNameAndDebug(patchSwitch, 0);
    addToRunTimeSelectionTable(fvMeshGIBChanger, patchSwitch, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshGIBChangers::patchSwitch::patchSwitch(fvMesh& mesh)
:
    fvMeshGIBChanger(mesh),
    dynamicMeshCoeffs_(dict().optionalSubDict(typeName + "Coeffs")),
    pSwitchPtr_()
{
    pSwitchPtr_.setSize(dynamicMeshCoeffs_.size());

    label zoneI = 0;

    forAllConstIter(dictionary, dynamicMeshCoeffs_, iter)
    {
        if (iter().isDict())
        {
            const dictionary& subDict = iter().dict();

            pSwitchPtr_.set(zoneI, GIBSwitch::New(mesh, subDict));

            zoneI++;
        }
    }
}


Foam::fvMeshGIBChangers::patchSwitch::patchSwitch
(
    fvMesh& mesh,
    const dictionary dict
)
:
    fvMeshGIBChanger(mesh),
    dynamicMeshCoeffs_(dict.optionalSubDict(typeName + "Coeffs")),
    pSwitchPtr_()
{
    pSwitchPtr_.setSize(dynamicMeshCoeffs_.size());

    label zoneI = 0;

    forAllConstIter(dictionary, dynamicMeshCoeffs_, iter)
    {
        if (iter().isDict())
        {
            const dictionary& subDict = iter().dict();

            pSwitchPtr_.set(zoneI, GIBSwitch::New(mesh, subDict));

            zoneI++;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshGIBChangers::patchSwitch::update()
{
    bool meshChanged = false;

    forAll(pSwitchPtr_, psI)
    {
        bool changedI = pSwitchPtr_[psI].update();
        if (changedI == true)
        {
            meshChanged = true;
        }
    }

    return meshChanged;
}


// ************************************************************************* //
