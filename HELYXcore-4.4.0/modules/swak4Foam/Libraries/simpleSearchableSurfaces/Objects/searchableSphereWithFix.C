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
    (c) 1991-2008 OpenCFD Ltd.

Contributors/Copyright:
    2009, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#include "searchableSphereWithFix.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(searchableSphereWithFix, 0);
addToRunTimeSelectionTable(searchableSurface, searchableSphereWithFix, dict);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSphereWithFix::searchableSphereWithFix
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSphere(
        io,
        dict
    ),
    centre_(dict.lookup("centre")),
    radius_(dict.lookup<scalar>("radius"))
{
    WarningIn("searchableSphereWithFix::searchableSphereWithFix")
        << "This class is only a workaround for the original"
            << " searchableSphere-implementation which cause a division "
            << "by zero and is present in OF up to 2.2. Should not be used "
            << "for versions where this is fixed"
            << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableSphereWithFix::~searchableSphereWithFix()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::pointIndexHit Foam::searchableSphereWithFix::findNearest
(
    const point& sample,
    const scalar nearestDistSqr
) const
{
    pointIndexHit info(false, sample, -1);

    const vector n(sample-centre_);
    scalar magN = mag(n);

    if (nearestDistSqr > sqr(magN-radius_))
    {
        if (magN < ROOTVSMALL)
        {
            //            info.rawPoint() = centre_ + vector(1,0,0)/magN*radius_;
            info.rawPoint() = centre_ + vector(1,0,0)/ROOTVSMALL*radius_;
        }
        else
        {
            info.rawPoint() = centre_ + n/magN*radius_;
        }
        info.setHit();
        info.setIndex(0);
    }

    return info;
}

void Foam::searchableSphereWithFix::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info,
    const bool threaded /*= false*/
) const
{
    info.setSize(samples.size());

    forAll(samples, i)
    {
        info[i] = findNearest(samples[i], nearestDistSqr[i]);
    }
}

// ************************************************************************* //
