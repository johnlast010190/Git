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
    fieldDistributionFunctionObject

Description
    Distribution of a field on the interior

SourceFiles
    fieldDistributionFunctionObject.C

Contributors/Copyright:
    2008-2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef fieldDistributionFunctionObject_H
#define fieldDistributionFunctionObject_H

#include "distributionFunctionObject/distributionFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class fieldDistributionFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class fieldDistributionFunctionObject
:
    public distributionFunctionObject
{
    // Private Data

        word fieldName_;


    // Private Member Functions

        template <typename T>
        void getDistributionInternal(autoPtr<SimpleDistribution<T>>& dist);


protected:

    // Protected Member Functions

        word dirName();

        virtual word baseName();

        virtual void getDistribution();

        const word& fieldName() const
        {
            return fieldName_;
        }


public:

    //- Runtime type information
    TypeName("fieldDistribution");


    // Constructors

        //- Construct from components
        fieldDistributionFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        fieldDistributionFunctionObject
        (
            const fieldDistributionFunctionObject&
        ) = delete;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fieldDistributionFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "fieldDistributionFunctionObjectTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
