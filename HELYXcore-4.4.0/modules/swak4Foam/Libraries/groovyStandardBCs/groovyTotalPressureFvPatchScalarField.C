/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) ICE Stroemungsfoschungs GmbH
    (c) 1991-2010 OpenCFD Ltd.

Contributors/Copyright:
    2011, 2013-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "groovyTotalPressureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::groovyTotalPressureFvPatchScalarField::groovyTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    totalPressureFvPatchScalarField(p, iF),
    p0Expression_("0"),
    driver_(this->patch())
{}


Foam::groovyTotalPressureFvPatchScalarField::groovyTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    totalPressureFvPatchScalarField(p, iF,dict),
    p0Expression_(
        dict.lookup("p0Expression"),
        dict
    ),
    driver_(dict,this->patch(),this->db())
{
}


Foam::groovyTotalPressureFvPatchScalarField::groovyTotalPressureFvPatchScalarField
(
    const groovyTotalPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    totalPressureFvPatchScalarField(ptf, p, iF, mapper),
    p0Expression_(ptf.p0Expression_),
    driver_(this->patch(),ptf.driver_)
{}


Foam::groovyTotalPressureFvPatchScalarField::groovyTotalPressureFvPatchScalarField
(
    const groovyTotalPressureFvPatchScalarField& tppsf
)
:
    totalPressureFvPatchScalarField(tppsf),
    p0Expression_(tppsf.p0Expression_),
    driver_(this->patch(),tppsf.driver_)
{}


Foam::groovyTotalPressureFvPatchScalarField::groovyTotalPressureFvPatchScalarField
(
    const groovyTotalPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    totalPressureFvPatchScalarField(tppsf, iF),
    p0Expression_(tppsf.p0Expression_),
    driver_(this->patch(),tppsf.driver_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::groovyTotalPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    driver_.clearVariables();

    p0()=driver_.evaluate<scalar>(this->p0Expression_);

    totalPressureFvPatchScalarField::updateCoeffs();
}

void Foam::groovyTotalPressureFvPatchScalarField::write(Ostream& os) const
{
    totalPressureFvPatchScalarField::write(os);

    os.writeEntry("p0Expression", p0Expression_);
    driver_.writeCommon(os,debug);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        groovyTotalPressureFvPatchScalarField
    );
}

// ************************************************************************* //
