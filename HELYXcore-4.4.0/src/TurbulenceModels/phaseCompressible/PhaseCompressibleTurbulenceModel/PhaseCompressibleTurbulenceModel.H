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
    (c) 2013-2016 OpenFOAM Foundation

Class
    Foam::PhaseCompressibleTurbulenceModel

Description
    Templated abstract base class for multiphase compressible
    turbulence models.

SourceFiles
    PhaseCompressibleTurbulenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef PhaseCompressibleTurbulenceModel_H
#define PhaseCompressibleTurbulenceModel_H

#include "TurbulenceModel/TurbulenceModel.H"
#include "compressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class PhaseCompressibleTurbulenceModel Declaration
\*---------------------------------------------------------------------------*/

template<class TransportModel>
class PhaseCompressibleTurbulenceModel
:
    public TurbulenceModel
    <
        volScalarField,
        volScalarField,
        compressibleTurbulenceModel,
        TransportModel
    >
{

public:

    typedef volScalarField alphaField;
    typedef volScalarField rhoField;
    typedef TransportModel transportModel;


    // Constructors

        //- Construct
        PhaseCompressibleTurbulenceModel
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


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<PhaseCompressibleTurbulenceModel> New
        (
            const alphaField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transportModel,
            const word& propertiesName = turbulenceModel::propertiesName
        );


    //- Destructor
    virtual ~PhaseCompressibleTurbulenceModel()
    {}


    // Member Functions

        //- Return the laminar dynamic viscosity
        virtual const volScalarField& mu() const
        {
            return this->transport_.mu();
        }

        //- Return the turbulence dynamic viscosity
        virtual tmp<volScalarField> mut() const
        {
            return this->rho_*this->nut();
        }

        //- Return the turbulence dynamic viscosity on patch
        virtual tmp<scalarField> mut(const label patchi) const
        {
            return this->rho_.boundaryField()[patchi]*this->nut(patchi);
        }

        //- Return the effective dynamic viscosity
        virtual tmp<volScalarField> muEff() const
        {
            return mut() + mu();
        }

        //- Return the effective dynamic viscosity on patch
        virtual tmp<scalarField> muEff(const label patchi) const
        {
            return mut(patchi) + mu().boundaryField()[patchi];
        }

        //- Return the phase-pressure'
        // (derivative of phase-pressure w.r.t. phase-fraction)
        virtual tmp<volScalarField> pPrime() const;

        //- Return the face-phase-pressure'
        // (derivative of phase-pressure w.r.t. phase-fraction)
        virtual tmp<surfaceScalarField> pPrimef() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PhaseCompressibleTurbulenceModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
