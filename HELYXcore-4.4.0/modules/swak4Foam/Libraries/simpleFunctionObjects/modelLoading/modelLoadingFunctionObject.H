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
    modelLoadingFunctionObject

Description
    Basis for the other function objects the load and correct a model

SourceFiles
    modelLoadingFunctionObject.C

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef modelLoadingFunctionObject_H
#define modelLoadingFunctionObject_H

#include "simpleFunctionObject/simpleFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class modelLoadingFunctionObject Declaration
\*---------------------------------------------------------------------------*/

template <class ModelType>
class modelLoadingFunctionObject
:
    public simpleFunctionObject
{
    // Private Data

        //- The actual model we're working with
        autoPtr<ModelType> model_;

        //- Whether the model should be corrected
        bool correctModel_;

        //- Whether the model can be reloaded later
        bool allowReload_;

        //- Fail if a model of this type already exists
        bool failIfModelTypeExists_;


protected:

    // Protected Member Functions

        virtual autoPtr<ModelType> initModel() = 0;


public:

    // Constructors

        //- Construct from components
        modelLoadingFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        modelLoadingFunctionObject(const modelLoadingFunctionObject&) = delete;


    // Member Functions

        //- start is called at the start of the time-loop
        virtual inline bool start();

        //- Sample and write
        virtual void writeSimple();

#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
        bool read(const dictionary &);
#endif


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const modelLoadingFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#include "modelLoadingFunctionObject.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
