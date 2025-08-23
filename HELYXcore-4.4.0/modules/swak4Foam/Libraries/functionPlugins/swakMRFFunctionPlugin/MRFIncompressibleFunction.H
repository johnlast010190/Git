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
    Foam::MRFIncompressibleFunction

Description

SourceFiles
    MRFIncompressibleFunction.C

Contributors/Copyright:
    2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef MRFIncompressibleFunction_H
#define MRFIncompressibleFunction_H

#include "swakMRFPluginFunctionBasis.H"

#ifdef FOAM_HAS_IOMRFLIST

#include "plugins/FieldValuePluginFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class MRFIncompressibleFunction Declaration
\*---------------------------------------------------------------------------*/

template<class FieldType>
class MRFIncompressibleFunction
:
    public swakMRFPluginFunctionBasis
{
protected:

    // Protected Data

        autoPtr<FieldType> field_;


    // Protected Member Functions

        void doEvaluation();

        virtual void manipulate(FieldType &result) = 0;

        virtual void setArgument
        (
            label index,
            const string& content,
            const CommonValueExpressionDriver& driver
        );


public:

    // Constructors

        //- Construct from components
        MRFIncompressibleFunction
        (
            const FieldValueExpressionDriver& parentDriver,
            const word& name
        );

        //- Disallow default bitwise copy construction
        MRFIncompressibleFunction
        (
            const MRFIncompressibleFunction&
        ) = delete;


    //- Destructor
    virtual ~MRFIncompressibleFunction()
    {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const MRFIncompressibleFunction&) = delete;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // FOAM_HAS_IOMRFLIST

#endif

// ************************************************************************* //
