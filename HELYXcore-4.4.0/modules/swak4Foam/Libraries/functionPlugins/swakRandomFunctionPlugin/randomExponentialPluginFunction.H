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
    Foam::randomExponentialPluginFunction

Description

SourceFiles
    randomExponentialPluginFunction.C

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef randomExponentialPluginFunction_H
#define randomExponentialPluginFunction_H

#include "plugins/CommonPluginFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class randomExponentialPluginFunction Declaration
\*---------------------------------------------------------------------------*/

template <typename FType,typename DType>
class randomExponentialPluginFunction
:
    public FType
{
    // Private Member Functions

        //- Do the actual compuation
        void doEvaluationInternal(scalarField& f);


    // Private Data

        //- The halfLife of the exponential distribution
        scalar halfLife_;

        //- The seed for the random number generator
        label seed_;


public:

    //- Runtime type information
    TypeName("randomExponentialPluginFunction");


    // Constructors

        //- Construct from components
        randomExponentialPluginFunction
        (
            const DType& parentDriver,
            const word& name
        );

        randomExponentialPluginFunction
        (
            const randomExponentialPluginFunction&
        ) = delete;


    //- Destructor
    virtual ~randomExponentialPluginFunction()
    {}


    // Member Functions

        void doEvaluation();

        void setArgument(label index, const scalar&);

        void setArgument(label index, const label&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const randomExponentialPluginFunction&) = delete;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
