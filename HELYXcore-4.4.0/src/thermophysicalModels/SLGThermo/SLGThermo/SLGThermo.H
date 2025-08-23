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
    (c) 2011-2023 OpenFOAM Foundation

Class
    Foam::SLGThermo

Description
    Thermo package for (S)olids (L)iquids and (G)ases
    Takes reference to thermo package, and provides:
    - carrier : components of thermo - access to elemental properties
    - liquids : liquid components - access  to elemental properties
    - solids  : solid components - access  to elemental properties

    If thermo is not a multi-component thermo package, carrier is nullptr.
    Similarly, if no liquids or solids are specified, their respective
    pointers will also be nullptr.

    Registered to the mesh so that it can be looked-up

SourceFiles
    SLGThermo.C

\*---------------------------------------------------------------------------*/

#ifndef SLGThermo_H
#define SLGThermo_H

#include "db/regIOobject/regIOobject.H"
#include "fluidThermo/fluidThermo.H"
#include "mixtures/basicSpecieMixture/basicSpecieMixture.H"
#include "liquidProperties/liquidMixtureProperties/liquidMixtureProperties.H"
#include "solidProperties/solidMixtureProperties/solidMixtureProperties.H"
#include "materialModels/materialTables/materialTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class SLGThermo Declaration
\*---------------------------------------------------------------------------*/

class SLGThermo
:
    public regIOobject
{
    // Private data

        //- Thermo package
        fluidThermo& thermo_;

        //- Reference to the multi-component carrier phase thermo
        basicSpecieMixture* carrier_;

        //- Additional liquid properties data
        autoPtr<liquidMixtureProperties> liquids_;

        //- Additional solid properties data
        autoPtr<solidMixtureProperties> solids_;

        //- Pointer to the ha models
        matScalarUTable haModels_;

        //- Pointer to the hs models
        matScalarUTable hsModels_;

        //- Pointer to the hf models
        matScalarUTable hfModels_;

        //- Pointer to the Cp models
        matScalarUTable CpModels_;

        //- Pointer to the kappa models
        matScalarUTable kappaModels_;

        //- Pointer to the mu models
        matScalarUTable muModels_;


public:

    //- Runtime type information
    TypeName("SLGThermo");

    // Constructors

        //- Construct from mesh
        SLGThermo(const fvMesh& mesh, fluidThermo& thermo);


    //- Destructor
    virtual ~SLGThermo();


    // Member Functions

        //- Return reference to the thermo database
        const fluidThermo& thermo() const;

        scalar Cp(label speciei, scalar p, scalar T);

        scalar ha(label speciei, scalar p, scalar T);

        scalar kappa(label speciei, scalar p, scalar T);

        scalar hs(label speciei, scalar p, scalar T);

        scalar mu(label speciei, scalar p, scalar T);

        scalar hf(label speciei);

        //- Return reference to the gaseous components
        const basicSpecieMixture& carrier() const;

        //- Return reference to the global (additional) liquids
        const liquidMixtureProperties& liquids() const;

        //- Return reference to the global (additional) solids
        const solidMixtureProperties& solids() const;

        //- Index of carrier component
        label carrierId(const word& cmptName, bool allowNotFound = false) const;

        //- Index of liquid component
        label liquidId(const word& cmptName, bool allowNotFound = false) const;

        //- Index of solid component
        label solidId(const word& cmptName, bool allowNotFound = false) const;

        //- Thermo database has multi-component carrier flag
        bool hasMulticomponentCarrier() const;

        //- Thermo database has liquid components flag
        bool hasLiquids() const;

        //- Thermo database has solid components flag
        bool hasSolids() const;


        // IO

            bool writeData(Foam::Ostream&) const
            {
                return true;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
