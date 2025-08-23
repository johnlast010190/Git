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
    EvolveReactingCloudFunctionObject

Description
    Function object that creates and evolves a ReactingCloud at every timestep

SourceFiles
    EvolveReactingCloudFunctionObject.C

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#ifndef EvolveReactingCloudFunctionObject_H
#define EvolveReactingCloudFunctionObject_H

#include "EvolveCloudFunctionObject.H"

#include "include/swakCloudTypes.H"

#ifdef FOAM_REACTINGCLOUD_TEMPLATED
#include "BasicReactingCloud.H"
#else
#include "clouds/derived/basicReactingCloud/basicReactingCloud.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class EvolveReactingCloudFunctionObject Declaration
\*---------------------------------------------------------------------------*/

#ifdef FOAM_REACTINGCLOUD_TEMPLATED
template <class CloudType>
#endif
class EvolveReactingCloudFunctionObject
:
#ifdef FOAM_REACTINGCLOUD_TEMPLATED
    public EvolveCloudFunctionObject<CloudType>
#else
    public EvolveCloudFunctionObject<basicReactingCloud>
#endif
{
public:

    //- Runtime type information
    TypeName("EvolveReactingCloudFunctionObject");


    // Constructors

        //- Construct from components
        EvolveReactingCloudFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        EvolveReactingCloudFunctionObject
        (
            const EvolveReactingCloudFunctionObject&
        ) = delete;


    // Member Functions

        //- Start is called at the start of the time-loop
        virtual bool start();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const EvolveReactingCloudFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
