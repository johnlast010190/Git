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
    Foam::floodFillFromCellsPluginFunction

Description

SourceFiles
    floodFillFromCellsPluginFunction.C

Contributors/Copyright:
    2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef floodFillFromCellsPluginFunction_H
#define floodFillFromCellsPluginFunction_H

#include "floodFillGeneralPluginFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class floodFillFromCellsPluginFunction Declaration
\*---------------------------------------------------------------------------*/
class floodFillFromCellsPluginFunction
:
    public floodFillGeneralPluginFunction
{

    // Private Data

        autoPtr<volScalarField> blocked_;


protected:

    // Protected Member Data

        virtual void setArgument
        (
            label index,
            const string& content,
            const CommonValueExpressionDriver& driver
        );

        virtual void initFacesAndCells();


public:

    //- Runtime type information
    TypeName("floodFillFromCellsPluginFunction");


    // Constructors

        //- Construct from components
        floodFillFromCellsPluginFunction
        (
            const FieldValueExpressionDriver& parentDriver,
            const word& name
        );

        //- Disallow default bitwise copy construction
        floodFillFromCellsPluginFunction
        (
            const floodFillFromCellsPluginFunction&
        ) = delete;


    //- Destructor
    virtual ~floodFillFromCellsPluginFunction()
    {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const floodFillFromCellsPluginFunction&) = delete;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
