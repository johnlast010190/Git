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
    Foam::swakRadiationModelPluginFunction

Description

SourceFiles
    swakRadiationModelPluginFunction.C

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef swakRadiationModelPluginFunction_H
#define swakRadiationModelPluginFunction_H

#include "plugins/FieldValuePluginFunction.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "radiationModels/radiationModel/radiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class swakRadiationModelPluginFunction Declaration
\*---------------------------------------------------------------------------*/
class swakRadiationModelPluginFunction
:
    public FieldValuePluginFunction
{
protected:

    // Protected Member Functions

        static const radiationModel& radiationInternal
        (
            const fvMesh& reg
        );

        const radiationModel& radiation();


public:

    //- Runtime type information
    TypeName("swakRadiationModelPluginFunction");


    // Constructors

        //- Construct from components
        swakRadiationModelPluginFunction
        (
            const FieldValueExpressionDriver& parentDriver,
            const word& name,
            const word& returnValueType,
            const string& spec = ""
        );

        //- Disallow default bitwise copy construction
        swakRadiationModelPluginFunction
        (
            const swakRadiationModelPluginFunction&
        ) = delete;


    //- Destructor
    virtual ~swakRadiationModelPluginFunction() {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const swakRadiationModelPluginFunction&) = delete;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
