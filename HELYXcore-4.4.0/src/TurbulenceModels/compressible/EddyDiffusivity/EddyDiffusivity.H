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
    (c) 2015-2018 OpenFOAM Foundation
    (c) 2020 Engys Ltd.

Class
    Foam::EddyDiffusivity

Description
    Templated abstract base class for single-phase compressible
    turbulence models.

SourceFiles
    EddyDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef EddyDiffusivity_H
#define EddyDiffusivity_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class EddyDiffusivity Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class EddyDiffusivity
:
    public BasicTurbulenceModel
{

    // Private Member Functions

        //- helper function for alphat construction
        tmp<volScalarField> alphatReadIfPresent(const word& fieldName) const;

        //- helper function for alphat creation
        tmp<volScalarField> alphatAutoCreate() const;


protected:

    // Protected data

        // Model coefficients

            //- Turbulent Prandtl Number
            dimensionedScalar Prt_;

            //- Turbulent Schmidt Number
            dimensionedScalar Sct_;


        // Fields

            mutable autoPtr<volScalarField> alphat_;



    // Protected Member Functions

        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    // Constructors

        //- Construct
        EddyDiffusivity
        (
            const word& type,
            const alphaField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName
        );


    //- Destructor
    virtual ~EddyDiffusivity()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Turbulent thermal diffusivity for enthalpy [kg/m/s]
        virtual tmp<volScalarField> alphat() const;

        //- Turbulent thermal diffusivity for enthalpy for a patch [kg/m/s]
        virtual tmp<scalarField> alphat(const label patchi) const
        {
            return alphat()().boundaryField()[patchi];
        }

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture [J/m/s/K]
        virtual tmp<volScalarField> kappaEff() const
        {
            return this->transport_.kappaEff(alphat());
        }

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture for patch [J/m/s/K]
        virtual tmp<scalarField> kappaEff(const label patchi) const
        {
            return this->transport_.kappaEff(alphat(patchi), patchi);
        }

        //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
        virtual tmp<volScalarField> alphaEff() const
        {
            return this->transport_.alphaEff(alphat());
        }

        //- Effective thermal turbulent diffusivity of mixture
        //  for patch [kg/m/s]
        virtual tmp<scalarField> alphaEff(const label patchi) const
        {
            return this->transport_.alphaEff(alphat(patchi), patchi);
        }

        //- Return the turbulent mass diffusivity [kg/m/s]
        virtual tmp<volScalarField> Dmt() const;

        //- Return the turbulent mass diffusivity for a patch [kg/m/s]
        virtual tmp<scalarField> Dmt(const label patchi) const;

        //- Correct the turbulence thermal diffusivity for energy transport
        virtual void correctEnergyTransport();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "EddyDiffusivity/EddyDiffusivity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
