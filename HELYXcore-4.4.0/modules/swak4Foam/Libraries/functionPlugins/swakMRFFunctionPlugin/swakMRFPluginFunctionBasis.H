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
    Foam::swakMRFPluginFunction

Description
    Basis for all the MRF-functions. Basically makes sure that the MRF-data is found

SourceFiles
    swakMRFPluginFunction.C

Contributors/Copyright:
    2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef swakMRFPluginFunction_H
#define swakMRFPluginFunction_H

#include "include/swak.H"

#ifdef FOAM_HAS_IOMRFLIST

#include "plugins/FieldValuePluginFunction.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "sources/derived/MRFSource/IOMRFSourceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class swakMRFPluginFunction Declaration
\*---------------------------------------------------------------------------*/

class swakMRFPluginFunctionBasis
:
    public FieldValuePluginFunction
{
protected:

    // Protected Member Functions

        const fv::MRFSourceList& MRF();


public:

    // Constructors

        //- Construct from components
        swakMRFPluginFunctionBasis
        (
            const FieldValueExpressionDriver& parentDriver,
            const word& name,
            const word& returnValueType,
            const string& argumentSpecification
        );

        //- Disallow default bitwise copy construction
        swakMRFPluginFunctionBasis
        (
            const swakMRFPluginFunctionBasis&
        ) = delete;


    //- Destructor
    virtual ~swakMRFPluginFunctionBasis()
    {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const swakMRFPluginFunctionBasis&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // FOAM_HAS_IOMRFLIST

#endif

// ************************************************************************* //
