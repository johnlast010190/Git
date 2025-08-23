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
    (c) 2010-2024 Engys Ltd.
    (c) 2011-2019 OpenFOAM Foundation

Class
    Foam::RASModels::realizableKE

Group
    grpRASTurbulence

Description
    Realizable k-epsilon turbulence model for incompressible and compressible
    flows.

    References:
    \verbatim
        Shih, T. H., Liou, W. W., Shabbir, A., Yang, Z., & Zhu, J. (1994).
        A new k-epsilon eddy viscosity model for high Reynolds number
        turbulent flows: Model development and validation.
        NASA STI/Recon Technical Report N, 95, 11442.

        Shih, T. H., Liou, W. W., Shabbir, A., Yang, Z., & Zhu, J. (1995).
        A New k-epsilon Eddy Viscosity Model for High Reynolds Number
        Turbulent Flows.
        Computers and Fluids, 24(3), 227-238.
    \endverbatim

    The default model coefficients are
    \verbatim
        realizableKECoeffs
        {
            Cmu         0.09;
            A0          4.0;
            C2          1.9;
            sigmak      1.0;
            sigmaEps    1.2;
        }
    \endverbatim

SourceFiles
    realizableKE.C

\*---------------------------------------------------------------------------*/

#ifndef realizableKE_H
#define realizableKE_H

#include "RAS/RASModel/RASModel.H"
#include "eddyViscosity/eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                        Class realizableKE Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class realizableKE
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        using IOobject::modelName;

protected:

    // Protected data


        //  Segregated or Coupled
            Switch solveCoupled_;

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar A0_;
            dimensionedScalar C2_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;

            //- maximum allowed nut/nu
            scalar maxTurbViscRatio_;


        // Fields

            volScalarField k_;
            volScalarField epsilon_;


   // Protected Member Functions

        tmp<volScalarField> rCmu
        (
            const volTensorField& gradU,
            const volScalarField& S2,
            const volScalarField& magS
        );

        virtual void correctNut
        (
            const volTensorField& gradU,
            const volScalarField& S2,
            const volScalarField& magS
        );

        virtual void correctNut();
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;

        //- solve segregated
        virtual void segregatedCorrect();

        //- solve coupled
        virtual void coupledCorrect();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("realizableKE");

    // Constructors

        //- Construct from components
        realizableKE
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~realizableKE()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return volScalarField::New
            (
                "DkEff",
                (this->nut_/sigmak_ + this->nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return volScalarField::New
            (
                "DepsilonEff",
                (this->nut_/sigmaEps_ + this->nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return volScalarField::New("omega", epsilon_/(Cmu_*k_));
        }

        //- return turbulent length scale
        virtual tmp<volScalarField> L() const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "realizableKE.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
