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

#include "relaxationShapeCylindricalMoving.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "volMesh/volMesh.H"
#include "fields/GeometricFields/pointFields/pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationShapes
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(relaxationShapeCylindricalMoving, 0);
addToRunTimeSelectionTable
(
    relaxationShape,
    relaxationShapeCylindricalMoving,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


relaxationShapeCylindricalMoving::relaxationShapeCylindricalMoving
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    relaxationShape(subDictName, mesh_),
    coorFramePtr_
    (
        coeffDict_.found("referenceFrame")
      ? coordinateFrame::lookupNew(mesh_, coeffDict_)
      : nullptr
    ),    
    centre0_(vector(coeffDict_.lookup("centre"))),
    centre_(vector(coeffDict_.lookup("centre"))),
    rInner_(coeffDict_.lookup<scalar>("rInner")),
    rOuter_(coeffDict_.lookup<scalar>("rOuter")),
    domainCentre0_(gAverage(mesh_.C()))
{
    width_ = Foam::mag(rOuter_ - rInner_);

    if (mesh_.time().timeIndex() == 0)
    {
        centre0_ -= ( centre0_ & direction_ )*direction_;
    }
    else
    {
        if (coorFramePtr_)
        {
            Info<< "Restart from reference frames, no need to read points0 " << endl;
        }
        else
        {
            pointVectorField points0
            (
                IOobject
                (
                    "points0",
                    mesh_.time().findInstance(mesh_.meshDir(), "points0"),
                    polyMesh::meshSubDir,
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                pointMesh::New(mesh_)
            );
            centre0_ =
                - (centre0_ & direction_)*direction_ + gAverage(mesh_.points())
                - gAverage(points0.primitiveField());
        }
    }

    // Find computational cells inside the relaxation-shape
    findComputationalCells();

    // Computate the sigma coordinate
    computeSigmaCoordinate();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void relaxationShapeCylindricalMoving::findComputationalCells()
{
    const vectorField& cc = mesh_.C();

    cells_.setSize(5000);
    label count(0);

    // Translate the relaxation zone with the mesh
    if (coorFramePtr_)
    {
        const vector oldOrigin = coorFramePtr_->coorSys0().origin();
        const vector newOrigin = coorFramePtr_->coorSys().origin();

        centre_ = centre0_ + newOrigin - oldOrigin;
    }
    else
    {
        centre_ = centre0_ + gAverage(cc) - domainCentre0_;
    }

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


void relaxationShapeCylindricalMoving::computeSigmaCoordinate()
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


bool relaxationShapeCylindricalMoving::insideZone
(
    const label& celli
) const
{
    bool inside( false );

    vector cc( mesh_.C()[celli] );
    cc -= ( direction_ & cc )*direction_;

    scalar dist( Foam::mag( cc - centre_ ) );

    if (dist >= rInner_ && dist <= rOuter_)
    {
        inside = true;
    }

    return inside;
}


const pointField& relaxationShapeCylindricalMoving::pointSet()
{
    notImplemented("pointSet is not implemented for this shape");
}


scalar relaxationShapeCylindricalMoving::interpolation
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
