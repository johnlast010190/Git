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
    Foam::fv::reportAvailableFvOptionsFvOption

Description
    This fvOption does not contribute to the equation but only reports
    for *every* possible call that the call was done

SourceFiles
    reportAvailableFvOptionsFvOption.C

Contributors/Copyright:
    2014-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef reportAvailableFvOptionsFvOption_H
#define reportAvailableFvOptionsFvOption_H

#include "memory/autoPtr/autoPtr.H"
#include "sets/topoSetSource/topoSetSource.H"
#include "sets/topoSets/cellSet.H"
#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "cfdTools/general/fvOptions/fvOption.H"

#include "include/swakFvOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class reportAvailableFvOptionsFvOption Declaration
\*---------------------------------------------------------------------------*/

class reportAvailableFvOptionsFvOption
:
    public option
{

   // Private Member Functions

        void report(const string &message);


public:

    //- Runtime type information
    TypeName("reportAvailableFvOptions");


    // Constructors

        //- Construct from explicit source name and mesh
        reportAvailableFvOptionsFvOption
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        reportAvailableFvOptionsFvOption
        (
            const reportAvailableFvOptionsFvOption&
        ) = delete;


    //- Destructor
    ~reportAvailableFvOptionsFvOption();


    // Member Functions

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames);

        virtual bool alwaysApply() const;

        virtual label applyToField(const word &fieldName) const;

        //- Scalar
        virtual void correct(volScalarField& fld);

        //- Vector
        virtual void correct(volVectorField& fld);

        //- Spherical tensor
        virtual void correct(volSphericalTensorField& fld);

        //- Symmetric tensor
        virtual void correct(volSymmTensorField& fld);

        //- Tensor
        virtual void correct(volTensorField& fld);


        // Add explicit and implicit contributions to incompressible eqn

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
            virtual void setValue
            (
                fvMatrix<symmTensor>& eqn,
                const label fieldI
            );

            //- Tensor
            virtual void setValue(fvMatrix<tensor>& eqn, const label fieldI);


        // Constraints

            virtual void constrain
            (
                fvMatrix<scalar>& eqn,
                const label fieldi
            );

            virtual void constrain
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            virtual void constrain
            (
                fvMatrix<sphericalTensor>& eqn,
                const label fieldi
            );

            virtual void constrain
            (
                fvMatrix<symmTensor>& eqn,
                const label fieldi
            );

            virtual void constrain
            (
                fvMatrix<tensor>& eqn,
                const label fieldi
            );


        // Flux manipulations

            //- Make the given absolute flux relative
            virtual void makeRelative(surfaceScalarField& phi) const;

            //- Make the given absolute boundary flux relative
            virtual void makeRelative
            (
                FieldField<fvsPatchField, scalar>& phi
            ) const;

            //- Make the given absolute mass-flux relative
            virtual void makeRelative
            (
                const surfaceScalarField& rho,
                surfaceScalarField& phi
            ) const;

            //- Make the given relative flux absolute
            virtual void makeAbsolute(surfaceScalarField& phi) const;

            //- Make the given relative mass-flux absolute
            virtual void makeAbsolute
            (
                const surfaceScalarField& rho,
                surfaceScalarField& phi
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const reportAvailableFvOptionsFvOption&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
