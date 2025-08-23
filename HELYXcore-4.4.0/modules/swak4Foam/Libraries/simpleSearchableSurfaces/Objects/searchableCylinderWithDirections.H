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
    Foam::searchableCylinderWithDirections

Description
    Searching on sphere

SourceFiles
    searchableCylinderWithDirections.C

Contributors/Copyright:
    2009, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#ifndef searchableCylinderWithDirections_H
#define searchableCylinderWithDirections_H

#include "searchableSurfaces/searchableCylinder/searchableCylinder.H"

#include "include/swak.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class searchableCylinderWithDirections Declaration
\*---------------------------------------------------------------------------*/

class searchableCylinderWithDirections
:
    public searchableCylinder
{
    // Private Member Data

        //- Names of regions
        mutable wordList regions_;

        //- 'Left' point - duplicate from parent class as it is private
        const point point1_;

        //- 'Right' point - duplicate from parent class as it is private
        const point point2_;


    // Private Member Functions

        void recalcIndex(List<pointIndexHit>& info) const;


public:

    //- Runtime type information
    TypeName("searchableCylinderWithDirections");


    // Constructors

       //- Construct from dictionary (used by searchableSurface)
        searchableCylinderWithDirections
        (
            const IOobject& io,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        searchableCylinderWithDirections
        (
            const searchableCylinderWithDirections&
        ) = delete;


    //- Destructor
    virtual ~searchableCylinderWithDirections();


    // Member Functions

#ifdef FOAM_SEARCHABLE_SURF_USES_TMP
        virtual tmp<pointField>  coordinates() const;
#else
        virtual pointField coordinates() const;
#endif

        virtual const wordList& regions() const;

        virtual label size() const
        {
            return 3;
        }

        //- From a set of points and indices get the region
        virtual void getRegion
        (
            const List<pointIndexHit>&,
            labelList& region,
            const bool threaded = false
         ) const;


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

        bool writeData(Ostream&) const
        {
            notImplemented("searchableCylinderWithDirections::writeData(Ostream&) const");
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const searchableCylinderWithDirections&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
