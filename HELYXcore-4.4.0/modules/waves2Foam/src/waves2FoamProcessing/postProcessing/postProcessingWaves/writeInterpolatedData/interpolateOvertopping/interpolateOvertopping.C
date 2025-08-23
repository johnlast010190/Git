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

#include "interpolateOvertopping.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(interpolateOvertopping, 0);
addToRunTimeSelectionTable
(
    postProcessingWaves,
    interpolateOvertopping,
    postProcessingWaves
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


interpolateOvertopping::interpolateOvertopping
(
    const Time& rT,
    const dictionary& actionProp,
    const word& action
)
:
    postProcessingWaves( rT, actionProp, action )
{
}


interpolateOvertopping::~interpolateOvertopping()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void interpolateOvertopping::evaluate()
{
    List<std::pair<scalar, label>> timeLabel;
    wordList OTnames;
    List<vectorField> OTs;

    rawOvertopping rot( rT_, actionProperties_, actionType_ );

    rot.readOvertoppingData(timeLabel, OTnames, OTs);

    scalarField t = equidistantTime(timeLabel, actionProperties_);

    scalarField weights( t.size() );
    labelList leftData( t.size() );
    labelList rightData( t.size() );

    interpolationWeights(timeLabel, t, weights, leftData, rightData);

    vectorField output( weights.size(), vector::zero );

    forAll(OTs, OTI)
    {
        const vectorField& OT( OTs[OTI] );

        forAll(weights, ii)
        {
            output[ii] = weights[ii]*OT[leftData[ii]]
                + (1.0 - weights[ii] )*OT[rightData[ii]];
        }

        std::stringstream ss;
        ss << callName_ << "_" << OTI;

        writeIOVectorField( output, ss.str() );
    }

    std::stringstream ss;
    ss << callName_ << "_time";

    writeIOScalarField( t, ss.str() );

    writeNameDict( actionProperties_.lookup<scalar>("deltaT"), OTnames);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
