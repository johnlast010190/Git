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
    (c) 2015-2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::LESModels::WALE

Group
    grpLESTurbulence

Description
    The Wall-adapting local eddy-viscosity (WALE) SGS model.

    Reference:
    \verbatim
        Nicoud, F., & Ducros, F. (1999).
        Subgrid-scale stress modelling based on the square of the velocity
        gradient tensor.
        Flow, Turbulence and Combustion, 62(3), 183-200.
    \endverbatim

    The default model coefficients are
    \verbatim
        WALECoeffs
        {
            Ck                  0.094;
            Ce                  1.048;e
            Cw                  0.325;
        }
    \endverbatim

See also
    Foam::LESModels::Smagorinsky

SourceFiles
    WALE.C

\*---------------------------------------------------------------------------*/

#ifndef WALE_H
#define WALE_H

#include "LES/LESModel/LESModel.H"
#include "LES/LESeddyViscosity/LESeddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class WALE Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class WALE
:
    public LESeddyViscosity<BasicTurbulenceModel>
{
protected:

    // Protected data

        dimensionedScalar Cw_;


    // Protected Member Functions

        //- Return the deviatoric symmetric part of the square of the given
        //  velocity gradient field
        tmp<volSymmTensorField> Sd(const volTensorField& gradU) const;

        //- Return SGS kinetic energy
        //  calculated from the given velocity gradient
        tmp<volScalarField> k(const volTensorField& gradU) const;

        //- Update the SGS eddy-viscosity
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("WALE");


    // Constructors

        //- Construct from components
        WALE
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

        //- Disallow default bitwise copy construction
        WALE(const WALE&) = delete;


    //- Destructor
    virtual ~WALE()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k(fvc::grad(this->U_));
        }

        //- Correct Eddy-Viscosity and related properties
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const WALE&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "WALE.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
