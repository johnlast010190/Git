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
    Foam::mergeRegionsSearchableSurface

Description
    Makes the whole surface one region

SourceFiles
    mergeRegionsSearchableSurface.C

Contributors/Copyright:
    2009, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef mergeRegionsSearchableSurface_H
#define mergeRegionsSearchableSurface_H

#include "wrapperSearchableSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class mergeRegionsSearchableSurface Declaration
\*---------------------------------------------------------------------------*/

class mergeRegionsSearchableSurface
:
    public wrapperSearchableSurface
{
    // Private Data

        wordList regions_;
        pointField coordinates_;
        labelList toRegion_;


public:

    //- Runtime type information
    TypeName("mergeRegions");


    // Constructors

        //- Construct from dictionary (used by searchableSurface)
        mergeRegionsSearchableSurface
        (
            const IOobject& io,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        mergeRegionsSearchableSurface
        (
            const mergeRegionsSearchableSurface&
        ) = delete;


    //- Destructor
    virtual ~mergeRegionsSearchableSurface();


    // Member Functions

        virtual const wordList& regions() const;

        virtual label size() const { return regions().size(); }

#ifdef FOAM_SEARCHABLE_SURF_USES_TMP
        virtual tmp<pointField> coordinates() const;
#else
        virtual pointField coordinates() const;
#endif

        //- From a set of points and indices get the region
        virtual void getRegion
        (
            const List<pointIndexHit>&,
            labelList& region,
            const bool threaded = false
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const mergeRegionsSearchableSurface&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
