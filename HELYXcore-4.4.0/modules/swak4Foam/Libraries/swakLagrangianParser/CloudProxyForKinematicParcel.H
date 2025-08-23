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
    Foam::CloudProxyForKinematicParcel

Description

SourceFiles
    CloudProxyForKinematicParcel.C

Contributors/Copyright:
    2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef CloudProxyForKinematicParcel_H
#define CloudProxyForKinematicParcel_H

#include "CloudProxyForParticle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class CloudProxyForKinematicParcel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CloudProxyForKinematicParcel
:
    public CloudProxyForParticle<CloudType>
{
public:

    // Constructors

        //- Construct from components
        CloudProxyForKinematicParcel(const cloud&);

        //- Disallow default bitwise copy construction
        CloudProxyForKinematicParcel
        (
            const CloudProxyForKinematicParcel&
        ) = delete;


    //- Destructor
    virtual ~CloudProxyForKinematicParcel();


    // Member Functions

        tmp<scalarField> weights() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const CloudProxyForKinematicParcel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CloudProxyForKinematicParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
