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
    (c) 1991-2010 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::CloudProxyForReactingMultiphaseParcel

Description

SourceFiles
    CloudProxyForReactingMultiphaseParcel.C

Contributors/Copyright:
    2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef CloudProxyForReactingMultiphaseParcel_H
#define CloudProxyForReactingMultiphaseParcel_H

#include "CloudProxyForReactingParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class CloudProxyForReactingMultiphaseParcel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CloudProxyForReactingMultiphaseParcelNoComposition
:
    public CloudProxyForReactingParcelNoComposition<CloudType>
{
public:

    // Constructors

        //- Construct from components
        CloudProxyForReactingMultiphaseParcelNoComposition(const cloud&);

        //- Disallow default bitwise copy construction
        CloudProxyForReactingMultiphaseParcelNoComposition
        (
            const CloudProxyForReactingMultiphaseParcelNoComposition&
        ) = delete;


    //- Destructor
    virtual ~CloudProxyForReactingMultiphaseParcelNoComposition();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const CloudProxyForReactingMultiphaseParcelNoComposition&
        ) = delete;
};


template<class CloudType>
class CloudProxyForReactingMultiphaseParcel
:
    public CloudProxyForReactingMultiphaseParcelNoComposition<CloudType>
{
public:

    // Constructors

        //- Construct from components
        CloudProxyForReactingMultiphaseParcel(const cloud&);

        //- Disallow default bitwise copy construction
        CloudProxyForReactingMultiphaseParcel
        (
            const CloudProxyForReactingMultiphaseParcel&
        ) = delete;


    //- Destructor
    virtual ~CloudProxyForReactingMultiphaseParcel();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const CloudProxyForReactingMultiphaseParcel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CloudProxyForReactingMultiphaseParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
