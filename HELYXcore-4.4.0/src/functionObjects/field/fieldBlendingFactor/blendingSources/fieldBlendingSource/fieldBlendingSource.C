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
    (c) 2021 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "db/dictionary/dictionary.H"
#include "fieldBlendingSource.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
// Use surfaceInterpolate from finiteVolume, not this library!
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace blendingSources
{

defineTypeNameAndDebug(fieldBlendingSource, 0);
addToRunTimeSelectionTable
(blendingSource, fieldBlendingSource, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fieldBlendingSource::fieldBlendingSource
(
    const objectRegistry& obr,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    blendingSource(obr, mesh, dict),
    fieldName_(dict.lookup("fieldName"))
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<surfaceScalarField> fieldBlendingSource::sourceField()
{
    const volScalarField& field =
        obr_.lookupObject<volScalarField>(fieldName_);

    return fvc::interpolate(field);
}

// ************************************************************************* //

} // End namespace blendingSources
} // End namespace Foam

// ************************************************************************* //
