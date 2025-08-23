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
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::viscosityModels::BokilBewtra

Description
    An incompressible Bokil-Bewtra non-Newtonian viscosity model.

    The model is based on Bokil & Bewtra, 1972. The model blends linear
    to the fluid viscosity under a critical solid concentration value Xc:

    nu = nu_0 * 10^(aX) , X >= Xc
    nu = nu_0 * 10^(aX) * (Xc-X)/Xc + nu_l * (1 - (Xc-X)/Xc) , 0 < X < Xc

SourceFiles
    BokilBewtra.C

\*---------------------------------------------------------------------------*/

#ifndef BokilBewtra_H
#define BokilBewtra_H

#include "viscosityModels/viscosityModel/viscosityModel.H"
#include "dimensionedTypes/dimensionedScalar/dimensionedScalar.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class BokilBewtra Declaration
\*---------------------------------------------------------------------------*/

class BokilBewtra
:
    public viscosityModel
{
    // Private data

        dictionary BokilBewtraCoeffs_;

        // concentration field X
        word XName_;

        // mixture viscosity constant in [m^2/s]
        dimensionedScalar nu0_;

        // liquid viscosity in [m^2/s]
        dimensionedScalar nul_;

        // max viscosity in [m^2/s]
        dimensionedScalar numax_;

        // critical concentration in [g/L]
        dimensionedScalar Xc_;

        // exponent factor in [m^3/kg]
        dimensionedScalar a_;

        // used to convert X into [g/L]
        scalar fac_;

        volScalarField nu_;


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() const;


public:

    //- Runtime type information
    TypeName("BokilBewtra");


    // Constructors

        //- Construct from components
        BokilBewtra
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~BokilBewtra()
    {}


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Correct the laminar viscosity
        virtual void correct()
        {
            nu_ = calcNu();
        }

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
