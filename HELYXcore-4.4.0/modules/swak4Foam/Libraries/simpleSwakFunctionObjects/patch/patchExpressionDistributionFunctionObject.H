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

Class
    patchExpressionDistributionFunctionObject

Description
    Distribution of an expression on a patch

SourceFiles
    patchExpressionDistributionFunctionObject.C

Contributors/Copyright:
    2008-2011, 2013-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef patchExpressionDistributionFunctionObject_H
#define patchExpressionDistributionFunctionObject_H

#include "distributionFunctionObject/distributionFunctionObject.H"
#include "PatchValueExpressionDriver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class patchExpressionDistributionFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class patchExpressionDistributionFunctionObject
:
    public distributionFunctionObject
{
    // Private Data

        wordList patchNames_;
        labelList patchIDs_;

    // Private Member Functions

        template <typename T>
        void getDistributionInternal(autoPtr<SimpleDistribution<T>> &dist);


protected:

    // Protected Data

        //- The expression that is to be evaluated
        exprString expression_;

        //- The specification
        dictionary data_;

        PtrList<PatchValueExpressionDriver> drivers_;


    // Protected Member Functions

        word dirName();

        virtual word baseName();

        virtual void getDistribution();


public:

    //- Runtime type information
    TypeName("patchExpressionDistribution");


    // Constructors

        //- Construct from components
        patchExpressionDistributionFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        patchExpressionDistributionFunctionObject
        (
            const patchExpressionDistributionFunctionObject&
        ) = delete;


    // Member Functions

        bool start();

        void writeSimple();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const patchExpressionDistributionFunctionObject&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "patchExpressionDistributionFunctionObjectTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
