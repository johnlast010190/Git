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
    patchFieldFunctionObject

Description
    Basis for functionObjects that might work on every field on a patch

SourceFiles
    patchFieldFunctionObject.C

Contributors/Copyright:
    2008-2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef patchFieldFunctionObject_H
#define patchFieldFunctionObject_H

#include "patch/patchFunctionObject/patchFunctionObject.H"
#include "primitives/Tensor/tensor/tensor.H"
#include "primitives/SphericalTensor/sphericalTensor/sphericalTensor.H"
#include "primitives/SymmTensor/symmTensor/symmTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class patchFieldFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class patchFieldFunctionObject
:
    public patchFunctionObject
{
    // Private Data

        //- Names of fields to probe
        wordList fieldNames_;

        //- Indices (into fieldNames_) of scalar/vector/tensor fields
        wordList scalarFields_;
        wordList vectorFields_;
        wordList sphericalTensorFields_;
        wordList symmTensorFields_;
        wordList tensorFields_;


    // Private Member Functions

        //- Process and write field
        template<class T>
        void processAndWrite(const word& fieldName);

        //- Process and write field
        template<class T>
        void processAndWrite(const wordList& typeFields);

        template<class T>
        void findFields(wordList& typeFieldNames, boolList& foundFields);


protected:

    // Protected Member Functions

        //- Names of the files
        virtual wordList fileNames();

        //- This is clumsy, but virtual + templated
        //  member functions are not allowed
        virtual scalarField process(const word& fieldName, scalar preset) = 0;

        virtual Field<vector> process
        (
            const word& fieldName,
            vector preset
        ) = 0;

        virtual Field<sphericalTensor> process
        (
            const word& fieldName,
            sphericalTensor preset
        ) = 0;

        virtual Field<symmTensor> process
        (
            const word& fieldName,
            symmTensor preset
        ) = 0;

        virtual Field<tensor> process
        (
            const word& fieldName,
            tensor preset
        ) = 0;


public:

    //- Runtime type information
    TypeName("patchFieldFunctionObject");


    // Constructors

        //- Construct from components
        patchFieldFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        patchFieldFunctionObject(const patchFieldFunctionObject&) = delete;


    // Member Functions

        //- Start is called at the start of the time-loop
        virtual bool start();

        virtual void writeSimple();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const patchFieldFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "patchFieldFunctionObjectTemplates.C"
#endif


#endif

// ************************************************************************* //
