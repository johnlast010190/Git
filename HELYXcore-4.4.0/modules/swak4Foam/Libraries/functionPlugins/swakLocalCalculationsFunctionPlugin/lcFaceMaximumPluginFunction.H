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
    Foam::lcFaceMaximumPluginFunction

Description

SourceFiles
    lcFaceMaximumPluginFunction.C

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef lcFaceMaximumPluginFunction_H
#define lcFaceMaximumPluginFunction_H

#include "localCellCalculationFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class lcFaceMaximumPluginFunction Declaration
\*---------------------------------------------------------------------------*/
class lcFaceMaximumPluginFunction
:
    public localCellCalculationFunction
{
    // Private Data

        autoPtr<surfaceScalarField> original_;


protected:

    // Protected Member Functions

        virtual void doCellCalculation(volScalarField& f);

        virtual void setArgument
        (
            label index,
            const string& content,
            const CommonValueExpressionDriver& driver
        );


public:

    //- Runtime type information
    TypeName("lcFaceMaximumPluginFunction");


    // Constructors

        //- Construct from components
        lcFaceMaximumPluginFunction
        (
            const FieldValueExpressionDriver& parentDriver,
            const word& name
        );

        //- Disallow default bitwise copy construction
        lcFaceMaximumPluginFunction
        (
            const lcFaceMaximumPluginFunction&
        ) = delete;


    //- Destructor
    virtual ~lcFaceMaximumPluginFunction()
    {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const lcFaceMaximumPluginFunction&) = delete;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
