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
    (c) 2009 Icon CG Ltd.
    (c) 2010-2012 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/interpolatedCylindricalVelocity/interpolatedCylindricalVelocityFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/volFields/volFields.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "db/Time/Time.H"
#include "primitives/Scalar/scalar/scalar.H"


namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void interpolatedCylindricalVelocityFvPatchVectorField::mapTable()
{
    vectorField cf = patch().Cf();
    forAll(cf, fI)
    {
        vector axisToFace = cf[fI] - centre_;
        axisToFace -= (axisToFace & axis_) * axis_;

        scalar radius = mag(axisToFace);

        vector Uzrt = distribution_->value(radius);

        vector tangentialDir = (axis_ ^ axisToFace);
        tangentialDir /= stabilise(mag(tangentialDir), SMALL);

        vector Ugf = axis_ * Uzrt.x()
            + axisToFace/stabilise(mag(axisToFace), SMALL) * Uzrt.y()
            + tangentialDir * Uzrt.z();

        operator[](fI) = Ugf;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

interpolatedCylindricalVelocityFvPatchVectorField::
interpolatedCylindricalVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    distribution_(),
    phiName_("phi"),
    rhoName_("rho"),
    centre_(vector::zero),
    axis_(vector(1,0,0)),
    flowRate_()
{}


interpolatedCylindricalVelocityFvPatchVectorField::
interpolatedCylindricalVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF),
    distribution_(Function1<vector>::New("profile", dict)),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    centre_(dict.lookup("centre")),
    axis_(dict.lookup("axis")),
    flowRate_()
{
    axis_ /= mag(axis_);

    mapTable();

    if (dict.found("flowRate"))
    {
        flowRate_.reset(new scalar(dict.lookup<scalar>("flowRate")));
    }
}


interpolatedCylindricalVelocityFvPatchVectorField::
interpolatedCylindricalVelocityFvPatchVectorField
(
    const interpolatedCylindricalVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    distribution_(ptf.distribution_, false),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    centre_(ptf.centre_),
    axis_(ptf.axis_),
    flowRate_()
{
    if (ptf.flowRate_.valid())
    {
        flowRate_.reset(new scalar(ptf.flowRate_()));
    }
}


interpolatedCylindricalVelocityFvPatchVectorField::
interpolatedCylindricalVelocityFvPatchVectorField
(
    const interpolatedCylindricalVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    distribution_(ptf.distribution_, false),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    centre_(ptf.centre_),
    axis_(ptf.axis_),
    flowRate_()
{
    if (ptf.flowRate_.valid())
    {
        flowRate_.reset(new scalar(ptf.flowRate_()));
    }
}


interpolatedCylindricalVelocityFvPatchVectorField::
interpolatedCylindricalVelocityFvPatchVectorField
(
    const interpolatedCylindricalVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    distribution_(ptf.distribution_, false),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    centre_(ptf.centre_),
    axis_(ptf.axis_),
    flowRate_()
{
    if (ptf.flowRate_.valid())
    {
        flowRate_.reset(new scalar(ptf.flowRate_()));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void interpolatedCylindricalVelocityFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    mapTable();

    if (flowRate_.valid())
    {
        vectorField n(patch().nf());

        const surfaceScalarField& phi =
            db().lookupObject<surfaceScalarField>(phiName_);
        scalarField Uphi(*this & (-patch().Sf()));

        if (phi.dimensions() == dimVelocity*dimArea)
        {
            // volumetric flow-rate
            scalar currentFlowRate = gSum(Uphi);

            forceAssign(*this * flowRate_()/stabilise(currentFlowRate, SMALL));
        }
        else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
        {
            const fvPatchField<scalar>& rhop =
                patch().lookupPatchFieldInDb<volScalarField, scalar>(db(), rhoName_);

            // mass flow-rate
            scalar currentFlowRate = gSum(rhop*Uphi);

            forceAssign(*this * flowRate_()/stabilise(currentFlowRate, SMALL));
        }
        else
        {
            FatalErrorInFunction
                << "dimensions of " << phiName_ << " are incorrect" << nl
                << "    on patch " << this->patch().name()
                << " of field " << this->internalField().name()
                << " in file " << this->internalField().objectPath()
                << nl << exit(FatalError);
        }
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}

void interpolatedCylindricalVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchField<vector>::write(os);

    if (phiName_ != "phi")
    {
        os.writeEntry("phi", phiName_);
    }
    if (rhoName_ != "rho")
    {
        os.writeEntry("rho", rhoName_);
    }

    os.writeEntry("centre", centre_);
    os.writeEntry("axis", axis_);

    if (flowRate_.valid())
    {
        os.writeEntry("flowRate", flowRate_());
    }

    distribution_->writeData(os);

    this->writeEntry("value", os);
}

// ************************************************************************* //

makePatchTypeField
(
   fvPatchVectorField,
   interpolatedCylindricalVelocityFvPatchVectorField
);

} // End namespace Foam

// ************************************************************************* //
