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
    (c) 2022 Engys Ltd.
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::SolidThermo

Description
    Energy for a solid mixture

SourceFiles
    SolidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef SolidThermo_H
#define SolidThermo_H

#include "basicThermo/BasicThermo.H"
#include "materialModels/materialTables/materialTables.H"
#include "materialsSpecie/specieModels/specieModels.H"
#include "equationOfState/equationOfState/equationOfState.H"
#include "thermo/thermo/thermodynamics/thermodynamics.H"
#include "thermo/energyConversion/energyConversion/energyConversion.H"
#include "thermo/energy/energy/energy.H"
#include "thermo/energy/energyInverseTemperature/energyInverseTemperature.H"
#include "thermo/energyConversion/standardThermo/standardThermo.H"
#include "materialsSpecie/materialsSpecie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class SolidThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicSolidThermo, class MixtureType>
class SolidThermo
:
    public BasicThermo<BasicSolidThermo, MixtureType>
{
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();


    // Private Data

        //- Reference to coordinate frame
        mutable coordinateFrame* coorFramePtr_;


public:

    //- Runtime type information
    TypeName("SolidThermo");


    // Constructors

        //- Construct from obr and phase name
        SolidThermo
        (
            const objectRegistry&,
            const word& phaseName
        );

        //- Construct from obr, dictionary and phase name
        SolidThermo
        (
            const objectRegistry&,
            const dictionary&,
            const word& phaseName
        );

        //- Disallow default bitwise copy construction
        SolidThermo
        (
            const SolidThermo<BasicSolidThermo, MixtureType>&
        ) = delete;


    //- Destructor
    virtual ~SolidThermo();


    // Member functions

        //- Load coordinate system
        const coordinateSystem& loadCoor() const;

        //- Update properties
        virtual void correct();

        //- Anisotropic thermal conductivity [W/m/K]
        virtual tmp<volVectorField> Kappa() const;

        //- Anisotropic thermal conductivity [W/m/K]
        virtual tmp<vectorField> Kappa(const label patchi) const;

        //- Anisotropic thermal conductivity
        //  in the local coordinate system [W/m/K]
        virtual tmp<volSymmTensorField> KappaLocal() const;

        //- Anisotropic thermal conductivity for patch
        //  in the local coordinate system [W/m/K]
        virtual tmp<symmTensorField> KappaLocal(const label patchi) const;

        //- Anisotropic thermal diffusivity for enthalpy of mixture [kg/m/s]
        virtual tmp<volSymmTensorField> alphahLocal() const;

        //- Return true if thermal conductivity is isotropic
        virtual bool isotropic() const
        {
            return BasicThermo<BasicSolidThermo, MixtureType>::isotropic();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SolidThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
