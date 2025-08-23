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
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2023-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "motionSolvers/displacement/solidBody/pointPatchFields/derived/solidBodyMotionDisplacement/solidBodyMotionDisplacementPointPatchVectorField.H"
#include "fields/Fields/transformField/transformField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::solidBodyMotionDisplacementPointPatchVectorField::points0AndOffset()
{
    if (coorFramePtr_ && !coorFramePtr_->isDynamic())
    {
        coorFramePtr_->resetDynamic(true);
    }

    // Necessary to support old solvers
    // (helyxSolve has it's own updateStates)
    coordinateFrame::updateStates(patch().boundaryMesh().mesh().thisDb());

    const septernion transform =
        isSolidBody_
      ? SBMFPtr_().transformation()
      : coorFramePtr_->transformation();

    if (isSolidBody_)
    {
        DeprecationWarningInFunction
        (
            solidBodyMotionFunction::typeName,
            "entry for boundary",
            40200,
            "Please replace it by using referenceFrame."
        );
    }

    tmp<vectorField> p0AndOffset;

    const bool isIncremental =
        (coorFramePtr_ && coorFramePtr_->isIncrementalMotion())
     || (SBMFPtr_.valid() && SBMFPtr_().isIncrementalMotion());

    if (isIncremental)
    {
        if (!oldPoints_.valid())
        {
            oldPoints_.reset(new pointField(localPoints0()));
        }

        p0AndOffset = transformPoints(transform, oldPoints_());
        oldPoints_() = p0AndOffset.ref();
    }
    else
    {
        p0AndOffset = transformPoints(transform, localPoints0());
    }

    return p0AndOffset;
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionDisplacementPointPatchVectorField::
solidBodyMotionDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(p, iF),
    isSolidBody_(false),
    SBMFPtr_(),
    coorFramePtr_(nullptr)
{}


Foam::solidBodyMotionDisplacementPointPatchVectorField::
solidBodyMotionDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchVectorField(p, iF, dict, false),
    isSolidBody_(dict.found("solidBodyMotionFunction")),
    SBMFPtr_
    (
        isSolidBody_
      ? solidBodyMotionFunction::New(dict, this->db().time())
      : nullptr
    ),
    coorFramePtr_(nullptr)
{
    if (!isSolidBody_)
    {
        coorFramePtr_ =
            coordinateFrame::lookupNew
            (
                dynamic_cast<const fvMesh&>(this->db()),
                dict
            );
    }
    if (!dict.found("value"))
    {
        // Determine current local points and offset
        fixedValuePointPatchVectorField::forceAssign
        (
            points0AndOffset() - localPoints0()
        );
    }
}


Foam::solidBodyMotionDisplacementPointPatchVectorField::
solidBodyMotionDisplacementPointPatchVectorField
(
    const solidBodyMotionDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchVectorField(ptf, p, iF, mapper),
    isSolidBody_(ptf.isSolidBody_),
    SBMFPtr_(isSolidBody_ ? ptf.SBMFPtr_().clone().ptr() : nullptr),
    coorFramePtr_(!isSolidBody_ ? ptf.coorFramePtr_ : nullptr)
{
    fixedValuePointPatchVectorField::forceAssign
    (
        points0AndOffset() - localPoints0()
    );
}


Foam::solidBodyMotionDisplacementPointPatchVectorField::
solidBodyMotionDisplacementPointPatchVectorField
(
    const solidBodyMotionDisplacementPointPatchVectorField& ptf
)
:
    fixedValuePointPatchVectorField(ptf),
    isSolidBody_(ptf.isSolidBody_),
    SBMFPtr_(isSolidBody_ ? ptf.SBMFPtr_().clone().ptr() : nullptr),
    coorFramePtr_(!isSolidBody_ ? ptf.coorFramePtr_ : nullptr)
{
    if (coorFramePtr_)
    {
        coorFramePtr_->resetDynamic(true);
    }
}


Foam::solidBodyMotionDisplacementPointPatchVectorField::
solidBodyMotionDisplacementPointPatchVectorField
(
    const solidBodyMotionDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(ptf, iF),
    isSolidBody_(ptf.isSolidBody_),
    SBMFPtr_(isSolidBody_ ? ptf.SBMFPtr_().clone().ptr() : nullptr),
    coorFramePtr_(!isSolidBody_ ? ptf.coorFramePtr_ : nullptr)
{
    fixedValuePointPatchVectorField::forceAssign
    (
        points0AndOffset() - localPoints0()
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::pointField&
Foam::solidBodyMotionDisplacementPointPatchVectorField::localPoints0() const
{
    if (!localPoints0Ptr_.valid())
    {
        pointIOField points0
        (
            IOobject
            (
                "points",
                this->db().time().constant(),
                polyMesh::meshSubDir,
                this->db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        localPoints0Ptr_.reset(new pointField(points0, patch().meshPoints()));
    }
    return localPoints0Ptr_();
}


void Foam::solidBodyMotionDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Determine current local points and offset
    fixedValuePointPatchVectorField::forceAssign
    (
        points0AndOffset() - localPoints0()
    );

    fixedValuePointPatchVectorField::updateCoeffs();
}


void Foam::solidBodyMotionDisplacementPointPatchVectorField::write
(
    Ostream& os
) const
{
    // Note: write value
    fixedValuePointPatchVectorField::write(os);

    if (isSolidBody_)
    {
        os.writeEntry(solidBodyMotionFunction::typeName, SBMFPtr_->type());
        os  << indent << word(SBMFPtr_->type() + "Coeffs");
        SBMFPtr_->writeData(os);
    }
    else
    {
        os.writeEntry("referenceFrame", coorFramePtr_->name());
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePointPatchTypeField
    (
        pointPatchVectorField,
        solidBodyMotionDisplacementPointPatchVectorField
    );
}


// ************************************************************************* //
