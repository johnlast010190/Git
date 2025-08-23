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

#include "relaxationShapeSemiCylindrical.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationShapes
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(relaxationShapeSemiCylindrical, 0);
addToRunTimeSelectionTable
(
    relaxationShape,
    relaxationShapeSemiCylindrical,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


relaxationShapeSemiCylindrical::relaxationShapeSemiCylindrical
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    relaxationShape(subDictName, mesh_),

    centre_( vector(coeffDict_.lookup("centre")) ),
    zeroAngleDirection_( vector(coeffDict_.lookup("zeroAngleDirection")) ),
    rInner_( coeffDict_.lookup<scalar>("rInner") ),
    rOuter_( coeffDict_.lookup<scalar>("rOuter") ),
    angleStart_( coeffDict_.lookup<scalar>("angleStart") ),
    angleEnd_( coeffDict_.lookup<scalar>("angleEnd") )
{
    width_   = Foam::mag(rOuter_ - rInner_);
    centre_ -= ( centre_ & direction_ )*direction_;

    if (angleEnd_ > 180)
    {
        angleEnd_ -= 360;
    }
    if (angleEnd_ < -180)
    {
        angleEnd_ += 360;
    }
    if (angleStart_ > 180)
    {
        angleStart_ -= 360;
    }
    if (angleStart_ < -180)
    {
        angleStart_ += 360;
    }

    piHalfAngleDirection_ = ( zeroAngleDirection_ ^ direction_ );

    // Find computational cells inside the relaxation-shape
    findComputationalCells();

    // Computate the sigma coordinate
    computeSigmaCoordinate();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void relaxationShapeSemiCylindrical::findComputationalCells()
{
    const vectorField& cc = mesh_.C();

    cells_.setSize(5000);
    label count(0);

    forAll(cc, celli)
    {
        if (insideZone(celli))
        {
            cells_[count++] = celli;

            if (count == cells_.size())
            {
                cells_.setSize( static_cast<label>( count*1.1 ) );
            }
        }
    }

    cells_.setSize(count);
}


void relaxationShapeSemiCylindrical::computeSigmaCoordinate()
{
    const vectorField& C = mesh_.C();
    sigma_.setSize(cells_.size(), 0);

    forAll(cells_, celli)
    {
        vector cc( C[cells_[celli]] );
        cc -= ( ( cc & direction_ )*direction_ );

        sigma_[celli]  = ( Foam::mag( cc - centre_ ) - rInner_ )/width_;
    }
}


bool relaxationShapeSemiCylindrical::angleCheck
(
    const scalar& angle
) const
{
    if (angleStart_ < angleEnd_)
    {
        return ( angleStart_ <= angle && angle <= angleEnd_ );
    }
    else
    {
        return ( angleStart_ <= angle || angle <= angleEnd_ );
    }
}


bool relaxationShapeSemiCylindrical::insideZone
(
    const label& celli
) const
{
    bool inside( false );

    vector cc( mesh_.C()[celli] );
    cc -= ( direction_ & cc )*direction_;
    cc -= centre_;

    scalar dist( Foam::mag( cc ) );
    scalar angle = 180/PI_*Foam::atan2
        (
            cc & piHalfAngleDirection_,
            cc & zeroAngleDirection_
        );

    if (dist >= rInner_ && dist <= rOuter_ && angleCheck(angle))
    {
        inside = true;
    }

    return inside;
}


const pointField& relaxationShapeSemiCylindrical::pointSet()
{
    notImplemented("pointSet is not implemented for this shape");
}


scalar relaxationShapeSemiCylindrical::interpolation
(
    const scalarField& source,
    const point& p0
) const
{
    notImplemented("interpolation is not implemented for this shape");
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationShapes
} // End namespace Foam

// ************************************************************************* //
