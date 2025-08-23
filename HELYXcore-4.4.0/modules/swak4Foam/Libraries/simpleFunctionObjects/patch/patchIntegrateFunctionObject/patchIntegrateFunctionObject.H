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
    (c) 2024 Engys Ltd.

Class
    patchIntegrateFunctionObject

Description
    Integrates fields on the patches

SourceFiles
    patchIntegrateFunctionObject.C

Contributors/Copyright:
    2008-2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef patchIntegrateFunctionObject_H
#define patchIntegrateFunctionObject_H

#include "patch/patchFieldFunctionObject/patchFieldFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class patchIntegrateFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class patchIntegrateFunctionObject
:
    public patchFieldFunctionObject
{
    // Private Member Functions

        //- Integrate single field on all patches
        template <class T>
        Field<T> integrate(const word& fieldName, T preset) const;


protected:

    // Protected Member Functions

        word dirName();

        virtual scalarField process(const word& fieldName, scalar preset);

        virtual Field<vector> process(const word& fieldName, vector preset);

        virtual Field<sphericalTensor> process
        (
            const word& fieldName,
            sphericalTensor preset
        );

        virtual Field<symmTensor> process
        (
            const word& fieldName,
            symmTensor preset
        );

        virtual Field<tensor> process(const word& fieldName, tensor preset);


public:

    //- Runtime type information
    TypeName("patchIntegrate");


    // Constructors

        //- Construct from components
        patchIntegrateFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        patchIntegrateFunctionObject
        (
            const patchIntegrateFunctionObject&
        ) = delete;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const patchIntegrateFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "patchIntegrateFunctionObjectTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
