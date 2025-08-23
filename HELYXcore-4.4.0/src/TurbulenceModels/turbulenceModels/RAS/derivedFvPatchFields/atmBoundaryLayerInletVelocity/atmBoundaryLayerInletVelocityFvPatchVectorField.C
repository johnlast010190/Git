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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2022 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "atmBoundaryLayerInletVelocityFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

atmBoundaryLayerInletVelocityFvPatchVectorField::
atmBoundaryLayerInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    atmBoundaryLayer()
{}


atmBoundaryLayerInletVelocityFvPatchVectorField::
atmBoundaryLayerInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    atmBoundaryLayer(patch().Cf(), dict)
{
    vectorField::operator=(U(patch().Cf()));
}


atmBoundaryLayerInletVelocityFvPatchVectorField::
atmBoundaryLayerInletVelocityFvPatchVectorField
(
    const atmBoundaryLayerInletVelocityFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(pvf, p, iF, mapper),
    atmBoundaryLayer(pvf, mapper)
{}


atmBoundaryLayerInletVelocityFvPatchVectorField::
atmBoundaryLayerInletVelocityFvPatchVectorField
(
    const atmBoundaryLayerInletVelocityFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pvf, iF),
    atmBoundaryLayer(pvf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void atmBoundaryLayerInletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    atmBoundaryLayer::autoMap(m);
}


void atmBoundaryLayerInletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& pvf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(pvf, addr);

    const atmBoundaryLayerInletVelocityFvPatchVectorField& blpvf =
        refCast<const atmBoundaryLayerInletVelocityFvPatchVectorField>(pvf);

    atmBoundaryLayer::rmap(blpvf, addr);
}


void atmBoundaryLayerInletVelocityFvPatchVectorField::reset
(
    const fvPatchVectorField& pvf
)
{
    fixedValueFvPatchVectorField::reset(pvf);

    const atmBoundaryLayerInletVelocityFvPatchVectorField& blpvf =
        refCast<const atmBoundaryLayerInletVelocityFvPatchVectorField>(pvf);

    atmBoundaryLayer::reset(blpvf);
}


void atmBoundaryLayerInletVelocityFvPatchVectorField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchVectorField::autoMapGIB(mapper);
    atmBoundaryLayer::autoMapGIB(mapper);
}


void atmBoundaryLayerInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    atmBoundaryLayer::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    atmBoundaryLayerInletVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
