/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
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
    (c) 2012-2017 OpenFOAM Foundation
    (c) 2010-2024 Engys Ltd.

Class
    Foam::porosityModels::powerLaw

Description
    Power law porosity model, given by:

        \f[
            S = - \rho C_0 |U|^{(C_1 - 1)} U
        \f]

    where
    \vartable
        C_0      | model linear coefficient
        C_1      | model exponent coefficient
    \endvartable

SourceFiles
    powerLaw.C
    powerLawTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef powerLaw_H
#define powerLaw_H

#include "cfdTools/general/porosityModel/porosityModel/porosityModel.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace porosityModels
{

/*---------------------------------------------------------------------------*\
                      Class powerLaw Declaration
\*---------------------------------------------------------------------------*/

class powerLaw
:
    public porosityModel
{
protected:
    // Protected data

        //- C0 coefficient
        autoPtr<Function1<scalar>> C0_;

        //- C1 coefficient
        autoPtr<Function1<scalar>> C1_;

        //- Name of density field
        word rhoName_;


private:
    // Private Member Functions

        //- Apply resistance
        template<class RhoFieldType>
        void apply
        (
            scalarField& Udiag,
            vectorField& Usource,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;

        //- Apply resistance
        template<class RhoFieldType>
        void apply
        (
            tensorField& AU,
            vectorField& source,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;

        //- Apply resistance
        template<class RhoFieldType>
        void adjointApply
        (
            scalarField& Udiag,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& Uprimal
        ) const;


public:

    //- Runtime type information
    TypeName("powerLaw");


    // Constructors

        //- Constructor
        powerLaw
        (
            const word& name,
            const word& modelType,
            const objectRegistry& obr,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& cellZoneName
        );

        //- Disallow default bitwise copy construction
        powerLaw(const powerLaw&) = delete;

    //- Destructor
    virtual ~powerLaw();


    // Member Functions

        //- Transform the model data wrt mesh changes
        virtual void calcTransformModelData();

        //- Calculate the porosity force
        virtual void calcForce
        (
            const volVectorField& U,
            const volScalarField& rho,
            const volScalarField& mu,
            vectorField& force
        ) const;

        //- Add resistance
        virtual void correct(fvVectorMatrix& UEqn) const;

        virtual void correct(fvBlockMatrix<vector>& UEqn) const;

        //- Add resistance
        virtual void correct
        (
            fvVectorMatrix& UEqn,
            const volScalarField& rho,
            const volScalarField& mu
        ) const;

        virtual void correct
        (
            fvBlockMatrix<vector>& UEqn,
            const volScalarField& rho,
            const volScalarField& mu
        ) const;

        //- Add resistance
        virtual void correct
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU
        ) const;

        virtual void adjointCorrect
        (
            fvVectorMatrix& UaEqn,
            const volVectorField& Uprimal
        ) const;

        virtual void adjointCorrect
        (
            fvVectorMatrix& UaEqn,
            const volScalarField& rho,
            const volScalarField& mu,
            const volVectorField& Uprimal
        ) const;

        virtual void adjointCorrect
        (
            fvBlockMatrix<vector>&,
            const volVectorField&
        ) const;


        virtual void adjointCorrect
        (
            fvBlockMatrix<vector>&,
            const volScalarField&,
            const volScalarField&,
            const volVectorField&
        ) const
        {
            NotImplemented;
        }

    // I-O

        //- Write
        bool writeData(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const powerLaw&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace porosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cfdTools/general/porosityModel/powerLaw/powerLawTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
