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

#include "db/error/error.H"
#include "relaxationShape.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationShapes
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<relaxationShape> relaxationShape::New
(
    const word& subDictName,
    const fvMesh& mesh_
)
{
    word relaxationShapeTypeName;

    // Enclose the creation of the waveProp to ensure it is
    // deleted before the relaxationShape is created otherwise the dictionary
    // is entered in the database twice
    {
        const dictionary coeffDict_
        (
            (mesh_.thisDb().lookupObject<IOdictionary>("waveProperties"))
             .subDict(subDictName + "Coeffs")
             .subDict("relaxationZone")
        );

        relaxationShapeTypeName = coeffDict_.lookup<word>("relaxationShape");

    }

    const auto ctor =
        ctorTableLookup
        (
            "relaxation shape type",
            dictionaryConstructorTable_(),
            "relaxationShape"+relaxationShapeTypeName
        );
    return autoPtr<relaxationShape>(ctor(subDictName,mesh_));
}


autoPtr<relaxationShape> relaxationShape::New
(
    const word& subDictName,
    const word relaxationShapeTypeName,
    const fvMesh& mesh_
)
{
    const auto ctor =
        ctorTableLookup
        (
            "relaxation shape type",
            dictionaryConstructorTable_(),
            "relaxationShape"+relaxationShapeTypeName
        );
    return autoPtr<relaxationShape>(ctor(subDictName,mesh_));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationShapes
} // End namespace Foam

// ************************************************************************* //
