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
    (c) 2016-2024 Engys Ltd.
    (c) 2011-2012 OpenFOAM Foundation

Class
    Foam::LESModels::kOmegaSSTSASDES

Group
    grpLESTurbulence

Description
    kOmegaSSTSASDES LES turbulence model for incompressible flows
    based on:

    "Evaluation of the SST-SAS model: channel flow, asymmetric diffuser
    and axi-symmetric hill".
    European Conference on Computational Fluid Dynamics ECCOMAS CFD 2006.
    Lars Davidson


    The first term of the Qsas expression is corrected following:

    DESider A European Effort on Hybrid RANS-LES Modelling:
    Results of the European-Union Funded Project, 2004 - 2007
    (Notes on Numerical Fluid Mechanics and Multidisciplinary Design).
    Chapter 2, section 8 Formulation of the Scale-Adaptive Simulation (SAS)
    Model during the DESIDER Project. Published in Springer-Verlag Berlin
    Heidelberg 2009.
    F. R. Menter and Y. Egorov.


SourceFiles
    kOmegaSSTSASDES.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSSTSASDES_H
#define kOmegaSSTSASDES_H

#include "DES/kOmegaSSTDES/kOmegaSSTDES.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                        Class kOmegaSSTSASDES Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kOmegaSSTSASDES
:
    public kOmegaSSTDES<BasicTurbulenceModel>
{
protected:

    // Protected data

        // Model constants

            dimensionedScalar Cs_;
            dimensionedScalar alphaPhi_;
            dimensionedScalar zetaTilda2_;
            dimensionedScalar FSAS_;


    // Protected Member Functions

        tmp<volScalarField> Lvk2
        (
            const volScalarField& S2
        ) const;

        virtual tmp<fvScalarMatrix> Qsas
        (
            const volScalarField& S2,
            const volScalarField& gamma,
            const volScalarField& beta
        ) const;

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    //- Runtime type information
    TypeName("kOmegaSSTSASDES");


    // Constructors

        //- Construct from components
        kOmegaSSTSASDES
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
        kOmegaSSTSASDES(const kOmegaSSTSASDES&) = delete;


    //- Destructor
    virtual ~kOmegaSSTSASDES()
    {}


    // Member Functions

        //- Read LESProperties dictionary
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        kOmegaSSTSASDES& operator=(const kOmegaSSTSASDES&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#   include "kOmegaSSTSASDES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
