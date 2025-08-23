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
    (c) 2024 Engys Ltd.

Class
    Foam::rotateSearchableSurface

Description
    Searching on rotated surface

SourceFiles
    rotateSearchableSurface.C

Contributors/Copyright:
    2009, 2013-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef rotateSearchableSurface_H
#define rotateSearchableSurface_H

#include "transformationSearchableSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class rotateSearchableSurface Declaration
\*---------------------------------------------------------------------------*/

class rotateSearchableSurface
:
    public transformationSearchableSurface
{
    // Private Member Data

        tensor rotation_;
        tensor backRotation_;


protected:

    // Protected Member Functions

        //- Do the transformation for a point
        virtual point transform(const point &) const;

        //- Do the inverse transformation for a point
        virtual point inverseTransform(const point &) const;


public:

    //- Runtime type information
    TypeName("rotate");


    // Constructors

        //- Construct from components
        rotateSearchableSurface(const IOobject& io, const treeBoundBox& bb);

        //- Construct from dictionary (used by transformationSearchableSurface)
        rotateSearchableSurface
        (
            const IOobject& io,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        rotateSearchableSurface(const rotateSearchableSurface&) = delete;


    //- Destructor
    virtual ~rotateSearchableSurface();


    // Member Functions

        //- From a set of points and indices get the normal
        virtual void getNormal
        (
            const List<pointIndexHit>&,
            vectorField& normal
        ) const;

#ifdef FOAM_SEARCHABLE_SURF_NEEDS_BOUNDING_SPHERES
        virtual void boundingSpheres
        (
            pointField& centres,
            scalarField& radiusSqr
        ) const;
#endif


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const rotateSearchableSurface&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
