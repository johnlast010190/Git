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
    Foam::transformationSearchableSurface

Description
    Searching on transformed surface

SourceFiles
    transformationSearchableSurface.C

Contributors/Copyright:
    2009, 2013-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef transformationSearchableSurface_H
#define transformationSearchableSurface_H

#include "searchableSurfaces/searchableSurface/searchableSurface.H"

#include "include/swakSurfaceTypes.H"

#include "meshes/boundBox/boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class transformationSearchableSurface Declaration
\*---------------------------------------------------------------------------*/

class transformationSearchableSurface
:
    public searchableSurface
{
    // Private Member Data

        autoPtr<IOobject> rename_;

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

        //- Do the transformation for a point
        virtual point transform(const point&) const = 0;

        //- Do the inverse transformation for a point
        virtual point inverseTransform(const point&) const = 0;

        //- Transform a hit to the real world
        pointIndexHit transform(const pointIndexHit&) const;

        //- Transform a hit from the real world
        pointIndexHit inverseTransform(const pointIndexHit&) const;


public:

    //- Runtime type information
    TypeName("transformationSearchableSurface");


    // Constructors

        //- Construct from dictionary (used by searchableSurface)
        transformationSearchableSurface
        (
            const IOobject& io,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        transformationSearchableSurface
        (
            const transformationSearchableSurface&
        ) = delete;


    //- Destructor
    virtual ~transformationSearchableSurface();


    // Member Functions

        virtual label size() const
        {
            return delegate().size();
        }

#ifdef FOAM_SEARCHABLE_SURF_USES_TMP
        virtual tmp<pointField> coordinates() const;
#else
        virtual pointField coordinates() const;
#endif

#ifdef FOAM_SEARCHABLE_SURF_HAS_POINTS
        virtual tmp<pointField> points() const;
#endif

        virtual const wordList& regions() const;

        virtual bool overlaps(const boundBox &) const;

        //- Whether supports volume type below
        virtual bool hasVolumeType() const
        {
            return delegate().hasVolumeType();
        }


        // Multiple point queries.

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
                notImplemented("transformationSearchableSurface::writeData(Ostream&) const");
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const transformationSearchableSurface&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
