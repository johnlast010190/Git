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
    Foam::fv::reportMatrixFvOption

Description
    This fvOption does not contribute to the equation but only reports
    cumulative properties of the matrix

SourceFiles
    reportMatrixFvOption.C

Contributors/Copyright:
    2014-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef reportMatrixFvOption_H
#define reportMatrixFvOption_H

#include "memory/autoPtr/autoPtr.H"
#include "sets/topoSetSource/topoSetSource.H"
#include "sets/topoSets/cellSet.H"
#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "cfdTools/general/fvOptions/fvOption.H"
#include "db/functionObjects/functionObjectList/functionObjectList.H"
#include "include/swakFvOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                   Class reportMatrixFvOption Declaration
\*---------------------------------------------------------------------------*/

class reportMatrixFvOption
:
    public option
{
    // Private data

        bool doAddSup_;

        bool doSetValue_;


    // Private Member Functions

        template<class Type>
        void reportMatrix(const string& message, const fvMatrix<Type>& matrix);


public:

    //- Runtime type information
    TypeName("reportMatrix");


    // Constructors

        //- Construct from explicit source name and mesh
        reportMatrixFvOption
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        reportMatrixFvOption(const reportMatrixFvOption&) = delete;


    //- Destructor
    ~reportMatrixFvOption();


    // Member Functions

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames);

        //- Scalar
        virtual void addSup(fvMatrix<scalar>& eqn, const label fieldI);

        //- Vector
        virtual void addSup(fvMatrix<vector>& eqn, const label fieldI);

        //- Spherical tensor
        virtual void addSup(fvMatrix<symmTensor>& eqn, const label fieldI);

        //- Symmetric tensor
        virtual void addSup
        (
            fvMatrix<sphericalTensor>& eqn,
            const label fieldI
        );

        //- Tensor
        virtual void addSup(fvMatrix<tensor>& eqn, const label fieldI);


#ifdef FOAM_FVOPTION_HAS_ADDITIONAL_ADDSUP

    // Add explicit and implicit contributions to compressible equations

        //- Scalar
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldI
        );

        //- Vector
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldI
        );

        //- Spherical tensor
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<symmTensor>& eqn,
            const label fieldI
        );

        //- Symmetric tensor
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<sphericalTensor>& eqn,
            const label fieldI
        );

        //- Tensor
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<tensor>& eqn,
            const label fieldI
        );


    // Add explicit and implicit contributions to phase equations

        //- Scalar
        virtual void addSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldI
        );

        //- Vector
        virtual void addSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldI
        );

        //- Spherical tensor
        virtual void addSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvMatrix<symmTensor>& eqn,
            const label fieldI
        );

        //- Symmetric tensor
        virtual void addSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvMatrix<sphericalTensor>& eqn,
            const label fieldI
        );

        //- Tensor
        virtual void addSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvMatrix<tensor>& eqn,
            const label fieldI
        );

#endif

    // Set values directly

        //- Scalar
        virtual void setValue(fvMatrix<scalar>& eqn, const label fieldI);

        //- Vector
        virtual void setValue(fvMatrix<vector>& eqn, const label fieldI);

        //- Spherical tensor
        virtual void setValue
        (
            fvMatrix<sphericalTensor>& eqn,
            const label fieldI
        );

        //- Symmetric tensor
        virtual void setValue(fvMatrix<symmTensor>& eqn, const label fieldI);

        //- Tensor
        virtual void setValue(fvMatrix<tensor>& eqn, const label fieldI);


        //- in newer versions setValue was renamed
        virtual void constrain(fvMatrix<scalar>& eqn, const label fieldI)
        {
            setValue(eqn, fieldI);
        }

        virtual void constrain(fvMatrix<vector>& eqn, const label fieldI)
        {
            setValue(eqn, fieldI);
        }

        virtual void constrain(fvMatrix<tensor>& eqn, const label fieldI)
        {
            setValue(eqn, fieldI);
        }

        virtual void constrain(fvMatrix<symmTensor>& eqn, const label fieldI)
        {
            setValue(eqn, fieldI);
        }

        virtual void constrain
        (
            fvMatrix<sphericalTensor>& eqn,
            const label fieldI
        )
        {
            setValue(eqn, fieldI);
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const reportMatrixFvOption&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
