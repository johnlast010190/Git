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
    Foam::intersectSearchableSurface

Description
    Searching on rotated surface

SourceFiles
    intersectSearchableSurface.C

Contributors/Copyright:
    2009, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef intersectSearchableSurface_H
#define intersectSearchableSurface_H

#include "binaryOperationSearchableSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class intersectSearchableSurface Declaration
\*---------------------------------------------------------------------------*/

class intersectSearchableSurface
:
    public binaryOperationSearchableSurface
{
protected:

    // Protected Member Data

        virtual bool decidePoint
        (
            const hitWhom who,
            const bool inA,
            const bool inB
        ) const;

        //- Check whether the normal for a point
        //  on surface A should be reverted
        bool revertNormalA(const pointIndexHit&) const
        {
            return false;
        }

        //- Check whether the normal for a point
        //  on surface B should be reverted
        bool revertNormalB(const pointIndexHit&) const
        {
            return false;
        }


public:

    //- Runtime type information
    TypeName("intersect");


    // Constructors

        //- Construct from dictionary (used by transformationSearchableSurface)
        intersectSearchableSurface
        (
            const IOobject& io,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        intersectSearchableSurface(const intersectSearchableSurface&) = delete;


    //- Destructor
    virtual ~intersectSearchableSurface();


    // Member Functions

        void getVolumeType
        (
            const pointField& points,
            List<volumeType>& volType,
            const bool threaded = false
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const intersectSearchableSurface&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
