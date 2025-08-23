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

Class
    Foam::GeneralSurfacesPluginFunction

Description

SourceFiles
    GeneralSurfacesPluginFunction.C

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef GeneralSurfacesPluginFunction_H
#define GeneralSurfacesPluginFunction_H

#include "plugins/FieldValuePluginFunction.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "sampledSurface/sampledSurface/sampledSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class GeneralSurfacesPluginFunction Declaration
\*---------------------------------------------------------------------------*/
class GeneralSurfacesPluginFunction
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

        const sampledSurface& theSurface() const;

        const labelList& meshCells();


public:

    //- Runtime type information
    TypeName("GeneralSurfacesPluginFunction");


    // Constructors

        //- Construct from components
        GeneralSurfacesPluginFunction
        (
            const FieldValueExpressionDriver &parentDriver,
            const word &name,
            const word &resultType,
            const string &arguments
        );

        //- Disallow default bitwise copy construction
        GeneralSurfacesPluginFunction
        (
            const GeneralSurfacesPluginFunction&
        ) = delete;


    //- Destructor
    virtual ~GeneralSurfacesPluginFunction()
    {}


    // Member Functions

        //- Set the argument
        void setArgument
        (
            label index,
            const word &content
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const GeneralSurfacesPluginFunction&) = delete;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
