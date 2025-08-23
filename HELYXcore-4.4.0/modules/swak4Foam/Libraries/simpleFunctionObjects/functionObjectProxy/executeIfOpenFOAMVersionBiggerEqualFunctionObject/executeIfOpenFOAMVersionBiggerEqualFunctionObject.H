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
    executeIfOpenFOAMVersionBiggerEqualFunctionObject

Description
    This proxy-functionObject only executs if the OpenFOAM-version is bigger
    equal to the version specified by a triple
    - majorVersion
    - minorVersion
    - patchVersion

SourceFiles
    executeIfOpenFOAMVersionBiggerEqualFunctionObject.C

Contributors/Copyright:
    2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef executeIfOpenFOAMVersionBiggerEqualFunctionObject_H
#define executeIfOpenFOAMVersionBiggerEqualFunctionObject_H

#include "functionObjectProxy/conditionalFunctionObjectListProxy/conditionalFunctionObjectListProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class executeIfOpenFOAMVersionBiggerEqualFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class executeIfOpenFOAMVersionBiggerEqualFunctionObject
:
    public conditionalFunctionObjectListProxy
{
    // Private Data

        //- First part of the version
        label majorVersion_;

        //- Second part of the version
        label minorVersion_;

        //- Third (optional) part of the version
        label patchVersion_;


    // Private Member Functions

        //- Reread the data
        void readData(const dictionary& dict);

        //- Get a number as a label
        label toLabel(const string& v);


protected:

    //- the actual conition that determines whether the FO should be executed
    virtual bool condition();


public:

    //- Runtime type information
    TypeName("executeIfOpenFOAMVersionBiggerEqual");


    // Constructors

        //- Construct from components
        executeIfOpenFOAMVersionBiggerEqualFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        executeIfOpenFOAMVersionBiggerEqualFunctionObject
        (
            const executeIfOpenFOAMVersionBiggerEqualFunctionObject&
        ) = delete;


    // Member Functions

        bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const executeIfOpenFOAMVersionBiggerEqualFunctionObject&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
