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
    (c) 2024 Engys Ltd.

Class
    Foam::searchableSphereWithFix

Description
    Searching on sphere (fixes problem with the original implementation)

SourceFiles
    searchableSphereWithFix.C

Contributors/Copyright:
    2009, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#ifndef searchableSphereWithFix_H
#define searchableSphereWithFix_H

#include "searchableSurfaces/searchableSphere/searchableSphere.H"

#include "include/swak.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class searchableSphereWithFix Declaration
\*---------------------------------------------------------------------------*/

class searchableSphereWithFix
:
    public searchableSphere
{
    // Private Member Data

        //- Centre point
        const point centre_;

        //- Radius squared
        const scalar radius_;


    // Private Member Functions

        pointIndexHit findNearest
        (
            const point& sample,
            const scalar nearestDistSqr
        ) const;


public:

    //- Runtime type information
    TypeName("searchableSphereWithFix");


    // Constructors

       //- Construct from dictionary (used by searchableSurface)
        searchableSphereWithFix
        (
            const IOobject& io,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        searchableSphereWithFix(const searchableSphereWithFix&) = delete;


    //- Destructor
    virtual ~searchableSphereWithFix();


    // Member Functions

        virtual void findNearest
        (
            const pointField& sample,
            const scalarField& nearestDistSqr,
            List<pointIndexHit>&,
            const bool threaded = false
        ) const;


        bool writeData(Ostream&) const
        {
            notImplemented("searchableSphereWithFix::writeData(Ostream&) const");
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const searchableSphereWithFix&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
