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

#include "relaxationWeightCourantCorrection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationWeights
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(relaxationWeightCourantCorrection, 0);

// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //


scalar relaxationWeightCourantCorrection::courantNumber
(
    const surfaceScalarField& phi,
    const surfaceScalarField& magSf,
    const surfaceScalarField& deltaCoeffs,
    const cell& cc
) const
{
    scalar res(-1);
    scalar dt( mesh_.time().deltaT().value() );

    forAll(cc, facei)
    {
        if (cc[facei] < mesh_.nInternalFaces())
        {
            res = Foam::max
                (
                    res,
                    Foam::mag(phi[cc[facei]])/magSf[cc[facei]]
                   *deltaCoeffs[cc[facei]]*dt
                );
        }
    }

    return res;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


relaxationWeightCourantCorrection::relaxationWeightCourantCorrection
(
    const fvMesh& mesh,
    const dictionary& coeffDict
)
:
    mesh_(mesh),
    courantCorrection_
    (
        coeffDict.lookupOrDefault<Switch>("courantCorrection",false)
    )
{
}


relaxationWeightCourantCorrection::~relaxationWeightCourantCorrection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void relaxationWeightCourantCorrection::courantCorrection
(
    const labelList& cells,
    scalarField& weights
)
{
    if (courantCorrection_)
    {
        scalar maxCo( mesh_.time().controlDict().lookup<scalar>("maxCo") );

        const surfaceScalarField& phi
            (
                mesh_.thisDb().lookupObject<surfaceScalarField>("phi")
            );
        const surfaceScalarField& magSf( mesh_.magSf() );
        const surfaceScalarField& delta( mesh_.deltaCoeffs() );
        const cellList& cellFaces( mesh_.cells() );

        forAll(cells, celli)
        {
            label cellI( cells[celli] );

            scalar cn = courantNumber( phi, magSf, delta, cellFaces[cellI] );
            // Make a change in format to fit with the definition by
            // Sopheak Seng (SOSE)
            scalar w  = 1.0 - weights[celli];

            // The "1.0 - " from SOSE in front of Foam::pow is omitted to
            // correct for the local format change just above
            weights[celli] = Foam::pow( 1.0 - w, cn/maxCo );
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationWeights
} // End namespace Foam

// ************************************************************************* //
