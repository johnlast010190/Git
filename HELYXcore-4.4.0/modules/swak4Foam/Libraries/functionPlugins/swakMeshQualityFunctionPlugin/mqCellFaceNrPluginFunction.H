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
    Foam::mqCellFaceNrPluginFunction

Description

SourceFiles
    mqCellFaceNrPluginFunction.C

Contributors/Copyright:
    2012-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef mqCellFaceNrPluginFunction_H
#define mqCellFaceNrPluginFunction_H

#include "plugins/FieldValuePluginFunction.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class mqCellFaceNrPluginFunction Declaration
\*---------------------------------------------------------------------------*/
class mqCellFaceNrPluginFunction
:
    public FieldValuePluginFunction
{
public:

    //- Runtime type information
    TypeName("mqCellFaceNrPluginFunction");


    // Constructors

        //- Construct from components
        mqCellFaceNrPluginFunction
        (
            const FieldValueExpressionDriver& parentDriver,
            const word& name
        );

        //- Disallow default bitwise copy construction
        mqCellFaceNrPluginFunction
        (
            const mqCellFaceNrPluginFunction&
        ) = delete;


    //- Destructor
    virtual ~mqCellFaceNrPluginFunction()
    {}


    // Member Functions

        //- Perform the evaluation
        void doEvaluation();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const mqCellFaceNrPluginFunction&) = delete;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
