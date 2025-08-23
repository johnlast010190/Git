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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2010-2024 Engys Ltd.

Class
    Foam::LESModels::SpalartAllmarasDDES

Group
    grpDESTurbulence

Description
    SpalartAllmaras DDES turbulence model for incompressible and compressible
    flows

    Reference:
    \verbatim
        Spalart, P. R., Deck, S., Shur, M. L., Squires, K. D., Strelets, M. K.,
        & Travin, A. (2006).
        A new version of detached-eddy simulation, resistant to ambiguous grid
        densities.
        Theoretical and computational fluid dynamics, 20(3), 181-195.
    \endverbatim

SourceFiles
    SpalartAllmarasDDES.C

\*---------------------------------------------------------------------------*/

#ifndef SpalartAllmarasDDES_H
#define SpalartAllmarasDDES_H

#include "DES/SpalartAllmarasDES/SpalartAllmarasDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                     Class SpalartAllmarasDDES Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class SpalartAllmarasDDES
:
    public SpalartAllmarasDES<BasicTurbulenceModel>
{
    // Private Member Functions

        tmp<volScalarField> fd(const volScalarField& magGradU) const;

        tmp<volScalarField> rd(const volScalarField& magGradU) const;


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cd1_;
            dimensionedScalar Cd2_;


    // Protected Member Functions

        using IOobject::modelName;

        //- Length scale
        virtual tmp<volScalarField> dTilda
        (
            const volScalarField& chi,
            const volScalarField& fv1,
            const volTensorField& gradU
        ) const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("SpalartAllmarasDDES");


    // Constructors

        //- Construct from components
        SpalartAllmarasDDES
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
        SpalartAllmarasDDES(const SpalartAllmarasDDES&) = delete;


    //- Destructor
    virtual ~SpalartAllmarasDDES()
    {}


    // Member Functions

        //- Read from dictionary
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const SpalartAllmarasDDES&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SpalartAllmarasDDES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
