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
    (c) 2024 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "db/error/error.H"
#include "waveTheory.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


autoPtr<waveTheory> waveTheory::New
(
    const word& subDictName,
    const fvMesh& mesh_
)
{
    word waveTheoryTypeName;

    IOdictionary* waveProps =
        mesh_.thisDb().lookupObjectRefPtr<IOdictionary>("waveProperties");
    if (!waveProps)
    {
        waveProps =
            new IOdictionary
            (
                IOobject
                (
                    "waveProperties",
                    mesh_.time().constant(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );
        waveProps->store();
    }
    if (!mesh_.thisDb().foundObject<uniformDimensionedVectorField>("g"))
    {
        (
            new uniformDimensionedVectorField
            (
                IOobject
                (
                    "g",
                    mesh_.time().constant(),
                    mesh_,
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE
                )
            )
        )->store();
    }

    const dictionary coeffDict_(waveProps->subDict(subDictName + "Coeffs"));

    waveTheoryTypeName = coeffDict_.lookup<word>("waveType");

    const auto ctor =
        ctorTableLookup
        (
            "wave theory type",
            dictionaryConstructorTable_(),
            waveTheoryTypeName
        );
    return autoPtr<waveTheory>(ctor(subDictName, mesh_));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
