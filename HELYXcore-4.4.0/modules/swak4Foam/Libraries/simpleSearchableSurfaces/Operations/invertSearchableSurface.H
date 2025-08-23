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
    Foam::invertSearchableSurface

Description
    Invert the volume of the wrapped surface

SourceFiles
    invertSearchableSurface.C

Contributors/Copyright:
    2009, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef invertSearchableSurface_H
#define invertSearchableSurface_H

#include "wrapperSearchableSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class invertSearchableSurface Declaration
\*---------------------------------------------------------------------------*/

class invertSearchableSurface
:
    public wrapperSearchableSurface
{
public:

    //- Runtime type information
    TypeName("invert");


    // Constructors

        //- Construct from dictionary (used by searchableSurface)
        invertSearchableSurface
        (
            const IOobject& io,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        invertSearchableSurface(const invertSearchableSurface&) = delete;


    //- Destructor
    virtual ~invertSearchableSurface();


    // Member Functions

        void getVolumeType
        (
            const pointField& points,
            List<volumeType>& volType,
            const bool threaded = false
        ) const;

        //- From a set of points and indices get the normal
        void getNormal
        (
            const List<pointIndexHit>&,
            vectorField& normal
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const invertSearchableSurface&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
