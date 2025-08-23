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
    (c) 2023 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvsPatchFields/derived/polyFaces/polyFacesFvsPatchLabelField.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polyFacesFvsPatchLabelField::init()
{
    labelField::operator=(identity(patch().size()) + patch().start());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyFacesFvsPatchLabelField::polyFacesFvsPatchLabelField
(
    const fvPatch& p,
    const DimensionedField<label, surfaceMesh>& iF
)
:
    fvsPatchLabelField(p, iF)
{
    init();
}


Foam::polyFacesFvsPatchLabelField::polyFacesFvsPatchLabelField
(
    const fvPatch& p,
    const DimensionedField<label, surfaceMesh>& iF,
    const dictionary& dict
)
:
    fvsPatchLabelField(p, iF, dict, false)
{
    init();
}


Foam::polyFacesFvsPatchLabelField::polyFacesFvsPatchLabelField
(
    const polyFacesFvsPatchLabelField& ptf,
    const fvPatch& p,
    const DimensionedField<label, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvsPatchLabelField(ptf, p, iF, mapper, false)
{
    init();
}


Foam::polyFacesFvsPatchLabelField::polyFacesFvsPatchLabelField
(
    const polyFacesFvsPatchLabelField& ptf,
    const DimensionedField<label, surfaceMesh>& iF
)
:
    fvsPatchLabelField(ptf, iF)
{
    init();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeUntemplatedFvsPatchTypeField
    (
        fvsPatchLabelField,
        polyFacesFvsPatchLabelField
    );
}


// ************************************************************************* //
