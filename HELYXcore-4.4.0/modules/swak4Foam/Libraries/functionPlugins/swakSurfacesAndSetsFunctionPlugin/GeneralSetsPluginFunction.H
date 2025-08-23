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
    Foam::GeneralSetsPluginFunction

Description

SourceFiles
    GeneralSetsPluginFunction.C

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef GeneralSetsPluginFunction_H
#define GeneralSetsPluginFunction_H

#include "plugins/FieldValuePluginFunction.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "sampledSet/sampledSet/sampledSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class GeneralSetsPluginFunction Declaration
\*---------------------------------------------------------------------------*/
class GeneralSetsPluginFunction
:
    public FieldValuePluginFunction
{
    // Private Data

        word name_;
        autoPtr<labelList> meshCells_;


protected:

    // Protected Member Functions

        inline const word& name() const
        {
            return name_;
        }

        const sampledSet& theSet() const;

public:

    //- Runtime type information
    TypeName("GeneralSetsPluginFunction");


    // Constructors

        //- Construct from components
        GeneralSetsPluginFunction
        (
            const FieldValueExpressionDriver& parentDriver,
            const word& name,
            const word& resultType,
            const string& arguments
        );

        //- Disallow default bitwise copy construction
        GeneralSetsPluginFunction(const GeneralSetsPluginFunction&) = delete;


    //- Destructor
    virtual ~GeneralSetsPluginFunction()
    {}


    // Member Functions

        void setArgument(label index, const word &content);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const GeneralSetsPluginFunction&) = delete;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
