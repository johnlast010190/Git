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
    (c) 2022-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "referenceFrames/motionCoordinateFrame/motionCoordinateFrame.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fvMatrices/fvMatrices.H"
#include "solidBodyMotionFunctions/morphModesMotion/morphModesMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(motionCoordinateFrame, 0);
    addToRunTimeSelectionTable
    (
        coordinateFrame,
        motionCoordinateFrame,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::motionCoordinateFrame::dimensionCheck
(
    const dimensionSet& dim1,
    const dimensionSet& dim2
) const
{
    if (dim1 != dim2)
    {
        FatalErrorInFunction
            << "UEqn dimensions do not match " << dim2
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::motionCoordinateFrame::motionCoordinateFrame
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& frameName
)
:
    coordinateFrame(mesh, dict, frameName),
    motionFunction_(nullptr)
{
    const dictionary& motionDict = dict.subDict("motionFunction");
    motionFunction_ =
        solidBodyMotionFunction::New(mesh_, motionDict, frameName);

    // Should be defined from the motion so sub-dicts work correctly
    this->isIncrementalMotion() = motionFunction_().isIncrementalMotion();

    // If any of the motions needs to be incremental whole stuck will
    // be defined as incremental. This is inportant expecially for dynamic
    // nesting and non-motion frames that have dynamic parent frame.
    if (anyIncremental() || isIncrementalMotion())
    {
        isIncrementalMotion() = true;
        motionFunction_().isIncrementalMotion() = true;
        forAll(parents(), framei)
        {
            parents()[framei].isIncrementalMotion() = true;
            if (isA<motionCoordinateFrame>(parents()[framei]))
            {
                dynamic_cast<motionCoordinateFrame&>
                (
                    parents()[framei]
                ).motion().isIncrementalMotion() = true;
            }
        }
    }

    // Check if motion allows outer corrector motion
    this->outerCorrectorMotion() = motionFunction_().outerCorrectorMotion();

    // Initialize the coordinate frame state data
    // This dissalows the motion of the frame on construction but allows the
    // motion of morthModesMotion since it has to move into position.
    motionCoordinateFrame::updateState
    (
        !isA<solidBodyMotionFunctions::morphModesMotion>(motionFunction_())
    );
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::motionCoordinateFrame::~motionCoordinateFrame()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::motionCoordinateFrame::updateState(bool construction) const
{
    if (validParentFrame())
    {
        parentFrame().updateState(construction);
    }

    if (!isUpdated())
    {
        // First store old time data, then update
        storeOldTimes();

        // Update current coordinate system
        // If outerCorrectorMotion is true, then the motion function's
        // transformation function is incremental with respect to outer
        // correctors; however, rotation, velocity and acceleration are not
        const label index = nFrameCorrector();
        const label size = index + 1;
        transformations_.setSize(size);
        if (!outerCorrectorMotion() && index > 0)
        {
            transformations_[index] = septernion::I;
        }
        else
        {
            transformations_[index] = motionFunction_().transformation();
        }

        rotations_.setSize(size);
        rotations_[index] = motionFunction_().rotation();

        // Do not move the coordinate system if on construction of the frames
        // since this can cause issues with dynamic flag not being initialized
        // yet.
        if (!construction)
        {
            updateCoordinateSystem();
        }

        // Update the coordinate frame state data
        velocities_.setSize(size);
        velocities_[index] = motionFunction_().velocity();

        if (mesh_.thisDb().time().timeIndex() == 0)
        {
            const_cast<vectorTuple&>(oldTime().velocity()) =
                velocities_[index];
        }

        accelerations_.setSize(size);
        accelerations_[index] = motionFunction_().acceleration();

        updateIndex_ = obr_.time().timeIndex();

        if (obr_.time().writeTime())
        {
            coordinateFrameState::write();
        }
    }
}


Foam::tmp<Foam::volVectorField> Foam::motionCoordinateFrame::frameVelocity
(
    const volVectorField& positions,
    bool addParentFrames
) const
{
    const vector omega = Omega();
    const vector vel = velocity().first();

    // Frame solid body velocity velocity
    tmp<volVectorField> U
    (
        new volVectorField
        (
            IOobject
            (
                "Urf",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            (
                (
                    dimensionedVector("Omega", dimless/dimTime, omega)
                 ^ (
                        positions
                      - dimensionedVector("CofR", dimLength, CofR())
                    )
                )
              + dimensionedVector("Ulinear", dimVelocity, vel)
            )
        )
    );

    if (addParentFrame(addParentFrames))
    {
        U.ref() += parentFrame().frameVelocity(positions, addParentFrames);
    }

    return U;
}


Foam::tmp<Foam::vectorField> Foam::motionCoordinateFrame::frameVelocity
(
    const vectorField& positions,
    bool addParentFrames
) const
{
    tmp<vectorField> U((Omega() ^ (positions - CofR())) + velocity().first());

    if (addParentFrame(addParentFrames))
    {
        U.ref() += parentFrame().frameVelocity(positions, addParentFrames);
    }

    return U;
}


Foam::vector Foam::motionCoordinateFrame::frameVelocity
(
    const vector& position,
    bool addParentFrames
) const
{
    vector U((Omega() ^ (position - CofR())) + velocity().first());

    if (addParentFrame(addParentFrames))
    {
        U += parentFrame().frameVelocity(position, addParentFrames);
    }

    return U;
}


bool Foam::motionCoordinateFrame::writeData(Ostream& os) const
{
    os.writeEntry("origin", coorSys0().origin());
    os.writeEntry("axis", coorSys0().e3());
    return true;
}


// ************************************************************************* //
