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
    Foam::fvcSweepFunctionPlugin

Description

SourceFiles
    fvcSweepPluginfunction.C

Contributors/Copyright:
    2015-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef fvcSweepPluginfunction_H
#define fvcSweepPluginfunction_H

#include "plugins/FieldValuePluginFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class fvcSweepPluginFunction Declaration
\*---------------------------------------------------------------------------*/

class fvcSweepPluginFunction
:
    public FieldValuePluginFunction
{
// Private Data

        //- Number of layers
        label nLayers_;

        //- Difference of alpha
        scalar alphaDiff_;

        //- The field to be sweeped
        autoPtr<volScalarField> field_;

        //- The field to be sweeped with
        autoPtr<volScalarField> alpha_;


public:

    //- Runtime type information
    TypeName("fvcSweepPluginFunction");


    // Constructors

        //- Construct from components
        fvcSweepPluginFunction
        (
            const FieldValueExpressionDriver& parentDriver,
            const word& name
        );

        //- Disallow default bitwise copy construction
        fvcSweepPluginFunction(const fvcSweepPluginFunction&) = delete;


    //- Destructor
    virtual ~fvcSweepPluginFunction()
    {}


    // Member Functions

        void doEvaluation();

        void setArgument
        (
            label index,
            const string& content,
            const CommonValueExpressionDriver& driver
        );

        void setArgument(label index, const scalar&);

        void setArgument(label index, const label&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvcSweepPluginFunction&) = delete;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
