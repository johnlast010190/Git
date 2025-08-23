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

Contributors/Copyright:
    2009, 2013-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "rotateSearchableSurface.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "primitives/transform/transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(rotateSearchableSurface, 0);
addToRunTimeSelectionTable(searchableSurface, rotateSearchableSurface, dict);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rotateSearchableSurface::rotateSearchableSurface
(
    const IOobject& io,
    const dictionary& dict
)
:
    transformationSearchableSurface(io,dict)
{
    vector from(dict.lookup("rotateFrom"));
    vector to(dict.lookup("rotateTo"));

    if (mag(from)<SMALL || mag(to)<SMALL) {
        FatalErrorIn("rotateSearchableSurface::rotateSearchableSurface")
            << "Vector " << from << " or " << to << " close to zero"
                << endl
                << abort(FatalError);
    }

    from/=mag(from);
    to/=mag(to);

    rotation_ = rotationTensor(from,to);
    backRotation_ = rotationTensor(to,from);

#ifdef FOAM_SEARCHABLE_SURF_HAS_BOUND_METHOD
    pointField pts(delegate().bounds().points());
    forAll(pts,i) {
        pts[i]=transform(pts[i]);
    }
    bounds()=boundBox(pts);
#endif
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rotateSearchableSurface::~rotateSearchableSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::rotateSearchableSurface::transform(const point &p) const
{
    return Foam::transform(rotation_,p);
}

Foam::point Foam::rotateSearchableSurface::inverseTransform(const point &p) const
{
    return Foam::transform(backRotation_,p);
}

void Foam::rotateSearchableSurface::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    vectorField iNormal;

    transformationSearchableSurface::getNormal
        (
            info,
            iNormal
        );

    normal.setSize(iNormal.size());

    forAll(normal,i) {
        normal[i]=transform(iNormal[i]);
    }
}

#ifdef FOAM_SEARCHABLE_SURF_NEEDS_BOUNDING_SPHERES
void Foam::rotateSearchableSurface::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    delegate().boundingSpheres(
        centres,
        radiusSqr
    );
    forAll(centres,i) {
        centres[i]=Foam::transform(rotation_,centres[i]);
    }
}
#endif

// ************************************************************************* //
