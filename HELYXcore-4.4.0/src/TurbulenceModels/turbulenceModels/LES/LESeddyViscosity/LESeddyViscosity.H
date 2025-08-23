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
    (c) 2024 Engys Ltd.

Class
    Foam::LESModels::LESeddyViscosity

Group
    grpLESTurbulence

Description
    Eddy viscosity LES SGS model base class

SourceFiles
    LESeddyViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef LESeddyViscosity_H
#define LESeddyViscosity_H

#include "LES/LESModel/LESModel.H"
#include "eddyViscosity/eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                      Class LESeddyViscosity Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class LESeddyViscosity
:
    public eddyViscosity<LESModel<BasicTurbulenceModel>>
{
protected:

    // Protected data

        dimensionedScalar Ck_;
        dimensionedScalar Ce_;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    // Constructors

        //- Construct from components
        LESeddyViscosity
        (
            const word& type,
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName
        );

        //- Disallow default bitwise copy construction
        LESeddyViscosity(const LESeddyViscosity&) = delete;


    //- Destructor
    virtual ~LESeddyViscosity()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return sub-grid disipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Return sub-grid specific dissipation rate
        virtual tmp<volScalarField> omega() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LESeddyViscosity&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LESeddyViscosity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
