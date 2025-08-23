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
    Foam::wrapperSearchableSurface

Description
    Searching on transformed surface

SourceFiles
    wrapperSearchableSurface.C

Contributors/Copyright:
    2009, 2013-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef wrapperSearchableSurface_H
#define wrapperSearchableSurface_H

#include "searchableSurfaces/searchableSurface/searchableSurface.H"

#include "include/swakSurfaceTypes.H"

#include "meshes/boundBox/boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class wrapperSearchableSurface Declaration
\*---------------------------------------------------------------------------*/

class wrapperSearchableSurface
:
    public searchableSurface
{
    // Private Member Data

        autoPtr<searchableSurface> delegate_;


protected:

    // Protected Member Functions

        searchableSurface& delegate()
        {
            return delegate_();
        }

        const searchableSurface& delegate() const
        {
            return delegate_();
        }


public:

    //- Runtime type information
    TypeName("wrapperSearchableSurface");


    // Constructors

        //- Construct from dictionary (used by searchableSurface)
        wrapperSearchableSurface
        (
            const IOobject& io,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        wrapperSearchableSurface(const wrapperSearchableSurface&) = delete;


    //- Destructor
    virtual ~wrapperSearchableSurface();


    // Member Functions

        virtual const wordList& regions() const;

        //- Whether supports volume type below
        virtual bool hasVolumeType() const
        {
            return delegate().hasVolumeType();
        }

        label size() const
        {
            return delegate().size();
        }


#ifdef FOAM_SEARCHABLE_SURF_USES_TMP
        virtual tmp<pointField> coordinates() const
#else
        virtual pointField coordinates() const
#endif
        {
            return delegate().coordinates();
        }

        // Multiple point queries

            virtual void findNearest
            (
                const pointField& sample,
                const scalarField& nearestDistSqr,
                List<pointIndexHit>&,
                const bool threaded = false
            ) const;

            virtual void findLine
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            virtual void findLineAny
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&,
                const bool threaded = false
            ) const;

            //- Get all intersections in order from start to end.
            virtual void findLineAll
            (
                const pointField& start,
                const pointField& end,
                List<List<pointIndexHit>>&
            ) const;

            //- From a set of points and indices get the region
            virtual void getRegion
            (
                const List<pointIndexHit>&,
                labelList& region,
                const bool threaded = false
            ) const;

            //- From a set of points and indices get the normal
            virtual void getNormal
            (
                const List<pointIndexHit>&,
                vectorField& normal
            ) const;

            void getVolumeType
            (
                const pointField& points,
                List<volumeType>& volType,
                const bool threaded = false
            ) const;


        // regIOobject implementation

            bool writeData(Ostream&) const
            {
                notImplemented("wrapperSearchableSurface::writeData(Ostream&) const");
            }

#ifdef FOAM_SEARCHABLE_SURF_HAS_POINTS
        virtual tmp<pointField> points() const
        {
            return delegate().points();
        }
#endif

        virtual bool overlaps(const boundBox&) const;

#ifdef FOAM_SEARCHABLE_SURF_NEEDS_BOUNDING_SPHERES
        virtual void boundingSpheres
        (
            pointField& centres,
            scalarField& radiusSqr
        ) const;
#endif


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const wrapperSearchableSurface&) = delete;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
