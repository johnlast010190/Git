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
    (c) 2017 OpenCFD Ltd.
    (c) 2011-2023 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd.

Class
    Foam::basicThermo

Description
    Abstract base-class for fluid and solid thermodynamic properties

SourceFiles
    basicThermo.C

\*---------------------------------------------------------------------------*/

#ifndef basicThermo_H
#define basicThermo_H

#include "fields/volFields/volFields.H"
#include "db/typeInfo/typeInfo.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "memory/autoPtr/autoPtr.H"
#include "primitives/strings/word/wordIOList.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class materialTables;

/*---------------------------------------------------------------------------*\
                        Class basicThermo Declaration
\*---------------------------------------------------------------------------*/
class basicThermo
{
protected:

    // Protected Member Functions

        //- Look up or construct field
        static volScalarField& lookupOrConstruct
        (
            const objectRegistry& obr,
            const word& name,
            const dimensionSet& dims = dimless,
            const scalar& defaultValue = 0
        );

        static volScalarField& lookupOrConstruct
        (
            const objectRegistry& obr,
            const char* name,
            const dimensionSet& dims = dimless,
            const scalar& defaultValue = 0
        );

        //- Return list of enthalpy/internal energy field boundary types
        //  by interrogating the temperature field boundary types
        wordList heBoundaryTypes();

        //- Return list of enthalpy/internal energy field boundary base types
        //  by interrogating the temperature field boundary types
        wordList heBoundaryBaseTypes();


public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;


    // Static Member data

        //- Dictionary name thermophysicalProperties/materialProperties
        static word dictName;

        //- Materials dict name
        static const word matDictName;


    // Static member functions

        //- Return the enthalpy/internal energy field boundary type
        //  by interrogating the temperature field boundary type
        static word heBoundaryType(const fvPatchScalarField& tpf);

        //- Return the enthalpy/internal energy field boundary base types
        //  by interrogating the temperature field boundary type
        static word heBoundaryBaseType(const fvPatchScalarField& tpf);

        //- Absolute temperature (T+TRef) [K] when using the thermo library;
        //  otherwise use looked-up T
        static tmp<volScalarField> TAbsIfFound
        (
            const objectRegistry& obr,
            const word& TName = "T"
        );

        //- Temperature offset (TRef) when using the thermo library,
        //  otherwise return 0
        static dimensionedScalar TRefIfFound(const objectRegistry& obr);


    //- Runtime type information
    TypeName("basicThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        basicThermo,
        objectRegistry,
        (const objectRegistry& obr, const word& phaseName),
        (obr, phaseName)
    );


    // Selectors

        //- Generic lookup for thermodynamics package thermoTypeName
        template<class Thermo, class Table>
        static typename Table::iterator lookupThermo
        (
            const dictionary& thermoTypeDict,
            Table& tablePtr,
            const int nCmpt,
            const char* cmptNames[],
            const word& thermoTypeName,
            const word& phaseName = word::null
        );

        //- Generic lookup for each of the related thermodynamics packages
        template<class Thermo, class Table>
        static typename Table::iterator lookupThermo
        (
            const dictionary& thermoDict,
            Table& tablePtr,
            const word& phaseName = word::null
        );

        //- Generic New for each of the related thermodynamics packages
        template<class Thermo>
        static autoPtr<Thermo> New
        (
            const objectRegistry&,
            const word& phaseName = word::null
        );

        //- Generic New for each of the related thermodynamics packages
        template<class Thermo>
        static autoPtr<Thermo> New
        (
            const objectRegistry&,
            const dictionary&,
            const word& phaseName = word::null
        );

        //- Specialisation of the Generic New for basicThermo
        static autoPtr<basicThermo> New
        (
            const objectRegistry&,
            const word& phaseName = word::null
        );

        //- Generic lookup for any thermo type
        static basicThermo* lookupPtr
        (
            const objectRegistry&,
            const word& phaseName = word::null
        );

        //- Generic lookup or create for any thermo type
        //  phaseName here is accepted only for thermophysicalProperties
        //  backward-compatilility. To look up a phasic material sub-model, use
        //  multiphaseThermo::lookupOrCreate(...) instead.
        static basicThermo& lookupOrCreate
        (
            const objectRegistry&,
            const word& phaseName = word::null
        );


    //- Destructor
    virtual ~basicThermo();


    // Member functions

        //- Mesh
        static const fvMesh& mesh(const objectRegistry& obr)
        {
            return fvSolutionRegistry::getMesh(obr);
        }

        //- Renowed access to object registry
        virtual const objectRegistry& db() const = 0;

        //- Loads either thermophysical or material properties dict
        static IOobject loadDictionary
        (
            const objectRegistry& obr,
            const word& phaseName,
            const bool isRegistred
        );

        //- Add phasic variable
        virtual void addPhasicVariable(const word& var) = 0;

        //- Is that variable phasic?
        virtual bool isPhasicVariable(const word& var) const = 0;

        //- Return phase name
        virtual const word& phaseName() const = 0;

        //- Name of a property for a given phase
        //  This function is overridden in the new material classes
        virtual word phasePropertyName
        (
            const word& name,
            const word& phaseName
        ) const
        {
            return name == "p" ? name : IOobject::groupName(name, phaseName);
        }

        //- Using phase name from thermo directly
        virtual word phasePropertyName(const word& name) const
        {
            return phasePropertyName(name, this->phaseName());
        }

        static const basicThermo& lookupThermo(const fvPatchScalarField& pf);

        //- Return dictionary for phase
        virtual const dictionary& phaseDict() const = 0;

        //- Properties dictionary
        virtual const IOdictionary& properties() const = 0;

        //- Non-const access the properties dictionary
        virtual IOdictionary& properties() = 0;

        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate(const string& app, const word&) const;

        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate(const string& app, const word&, const word&) const;

        //- Split name of thermo package into a list of the components names
        static wordList splitThermoName
        (
            const word& thermoName,
            const int nCmpt
        );

        //- Override to perform additional init after BasicThermo's init is run
        virtual void init() = 0;

        //- Update properties
        virtual void correct() = 0;

        //- Name of the thermo physics
        virtual word thermoName() const = 0;

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        virtual bool incompressible() const = 0;

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        virtual bool isochoric() const = 0;

        //- Return true if the Cp or Cv value (depending on selected energy
        //  type) is constant
        virtual bool isCpvConst() const = 0;

        //- It shouldn't be used from here. It is overriden in new material
        //  library, but we don't want the function to be forced in old thermo
        //  derived classes.
        virtual bool buoyant() const = 0;

        //- It the buoyant rho same as main rho?
        virtual bool distinctBuoyancy() const = 0;

        //- Should the dpdt term be included in the enthalpy equation
        virtual Switch dpdt() const = 0;

        //- Do we calculate T from he or vice-versa?
        virtual bool calculatesTFromhe() const
        {
            return true;
        }

        //- Set whether to solve for T from he or vice versa
        virtual void setToCalculateTFromhe(bool solveTFromhe) = 0;


        // Access to thermodynamic state variables

            //- Gauge pressure [Pa]
            //  Non-const access allowed for transport equations
            virtual volScalarField& p() = 0;

            //- Gauge pressure [Pa]
            virtual const volScalarField& p() const = 0;

            //- Return reference pressure for gauge pressure
            virtual dimensionedScalar pRef() const = 0;

            //- Return reference pressure for gauge pressure
            virtual scalar pRefValue() const = 0;

            //- Absolute pressure (p + pRef) [Pa]
            virtual tmp<volScalarField> pAbs() const = 0;

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const = 0;

            //- Density [kg/m^3] (density for use in buoyancy term)
            virtual tmp<volScalarField> buoyantRho() const = 0;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& he() = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& he() const = 0;

            //- Enthalpy/Internal energy/Temperature field
            //  (whichever is the input field to thermo)
            virtual const volScalarField& heT() const = 0;

            //- Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> he
            (
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const = 0;

            //- Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual tmp<volScalarField> hs() const = 0;

            //- Sensible enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> hs
            (
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Sensible enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> hs
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const = 0;

            //- Sensible enthalpy for patch [J/kg]
            virtual tmp<scalarField> hs
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Absolute enthalpy [J/kg]
            virtual tmp<volScalarField> ha() const = 0;

            //- Absolute enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> ha
            (
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Absolute enthalpy for patch [J/kg]
            virtual tmp<scalarField> ha
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Enthalpy of formation [J/kg]
            virtual tmp<volScalarField> hf() const = 0;

            //- Temperature from enthalpy/internal energy
            virtual tmp<volScalarField> The
            (
                const volScalarField& h,
                const volScalarField& p,
                const volScalarField& T0    // starting temperature
            ) const = 0;

            //- Temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> The
            (
                const scalarField& h,
                const scalarField& T0,      // starting temperature
                const label patchi
            ) const = 0;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual const volScalarField& Cp() const = 0;

            //- Heat capacity at constant volume [J/kg/K]
            virtual const volScalarField& Cv() const = 0;


        // Transport state

            //- Thermal conductivity of mixture [W/m/K]
            virtual const volScalarField& kappa() const = 0;


        // Derived thermodynamic properties

            //- Temperature
            virtual const volScalarField& T() const = 0;

            //- Return reference temperature to be added to T field
            virtual dimensionedScalar TRef() const = 0;

            //- Return reference temperature to be added to T field
            virtual scalar TRefValue() const = 0;

            //- Temperature relative to TRef
            //  Non-const access allowed for transport equations
            virtual volScalarField& T() = 0;

            //- Absolute temperature (T+TRef) [K]
            virtual tmp<volScalarField> TAbs() const = 0;


            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Heat capacity at constant pressure/volume [J/kg/K]
            virtual const volScalarField& Cpv() const = 0;

            //- Heat capacity at constant pressure/volume for patch [J/kg/K]
            virtual tmp<scalarField> Cpv
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Gamma = Cp/Cv []
            tmp<volScalarField> gamma() const;

            //- Gamma = Cp/Cv for patch []
            tmp<scalarField> gamma
            (
                const scalarField& T,
                const label patchi
            ) const;


    // Member operators

        //- Return material tables
        virtual materialTables& materials() const = 0;

        //- Return calculated model
        virtual tmp<volScalarField> operator[](const word& modelName) const = 0;
};


/*---------------------------------------------------------------------------*\
                 Class basicThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class basicThermo::implementation
:
    virtual public basicThermo,
    public IOdictionary
{
protected:

    // Protected Data

        //- Phase-name
        const word phaseName_;

        //- Material type
        const word materialType_;

        //- Reference pressure for gauge pressure
        scalar pRef_;

        //- Whether a reference pressure was provided
        bool pRefFound_;

        //- Reference temperature [K] to be added to temperature
        //  (typically used to allow T_ below to be in Celsius)
        scalar TRef_;

        //- Whether a reference temperature was provided
        bool TRefFound_;

        //- Temperature relative to TRef
        volScalarField& T_;

        //- Was the T field found rather than constructed?
        bool TFound_;

        //- Thermal conductivity [W/m/K]
        volScalarField kappa_;

        //- Should the dpdt term be included in the enthalpy equation
        Switch dpdt_;


    // Protected Member Functions

        //- Look up or construct field; if phaseName is valid, try to construct
        //  from the field without phase name if phasic one is not found
        volScalarField& lookupOrConstructPhasic
        (
            const objectRegistry& mesh,
            const word& name,
            bool& found,
            const dimensionSet& dims = dimless,
            const scalar& defaultValue = 0
        );

        //- Load pRef, etc from dict
        scalar readRefValue
        (
            const word& fieldName,
            const dimensionSet& dims,
            bool& refValueFound
        ) const;

        //- Lookup and check out field
        void lookupAndCheckout(const char* name) const;


public:

    // Static Member data

        //- Name of the thermophysical properties dictionary
        using basicThermo::dictName;

        //- Name of the materials dictionary
        using basicThermo::matDictName;


    // Constructors

        //- Construct from mesh and phase name
        implementation(const objectRegistry&, const word& phaseName);

        //- Construct from mesh, dictionary and phase name
        implementation
        (
            const objectRegistry&,
            const dictionary&,
            const word& phaseName
        );

        //- Disallow default bitwise copy construction
        implementation(const implementation&) = delete;



    //- Destructor
    virtual ~implementation();


    // Member Functions

        //- Override to perform additional init after BasicThermo's init is run
        virtual void init()
        {}

        //- Renowed access to object registry
        virtual const objectRegistry& db() const
        {
            return IOdictionary::db();
        }

        //- Add phasic variable
        virtual void addPhasicVariable(const word& var)
        {}

        //- Is that variable phasic?
        //  This function is overridden in the new material classes
        virtual bool isPhasicVariable(const word& var) const
        {
            return (var == "T");
        }

        //- Using phase name from thermo directly
        virtual word phasePropertyName
        (
            const word& name,
            const word& phaseName
        ) const
        {
            return basicThermo::phasePropertyName(name, phaseName);
        }

        //- Using phase name from thermo directly
        virtual word phasePropertyName(const word& name) const
        {
            return basicThermo::phasePropertyName(name);
        }

        //- Return dictionary for phase
        virtual const dictionary& phaseDict() const
        {
            if (phaseName_ == word::null)
            {
                return *this;
            }
            else
            {
                return this->subDict(phaseName_);
            }
        }

        //- Phase name
        const word& phaseName() const
        {
            return phaseName_;
        }

        //- Print fatal error for material properties requirement
        void materialPropertiesRequired() const
        {
            FatalErrorInFunction
                << "This solver requires the new material properties library."
                << nl << "Ensure 'system/materialProperties' is present."
                << endl << exit(FatalError);
        }

        //- Return true if the Cp or Cv value (depending on selected energy
        //  type) is constant
        virtual bool isCpvConst() const
        {
            materialPropertiesRequired();
            NotImplemented;
        }

        //- It shouldn't be used from here. It is overriden in new material
        //  library, but we don't want the function to be forced in old thermo
        //  derived classes.
        virtual bool buoyant() const
        {
            materialPropertiesRequired();
            NotImplemented;
        }

        //- It the buoyant rho same as main rho?
        virtual bool distinctBuoyancy() const
        {
            materialPropertiesRequired();
            NotImplemented;
        }

        //- Should the dpdt term be included in the enthalpy equation
        virtual Switch dpdt() const
        {
            return dpdt_;
        }

        //- Set whether to solve for T from he or vice versa
        virtual void setToCalculateTFromhe(bool solveTFromhe)
        {
            materialPropertiesRequired();
        }

        //- Return reference pressure for gauge pressure
        virtual dimensionedScalar pRef() const
        {
            return dimensionedScalar("pRef", dimPressure, pRef_);
        }

        //- Return reference pressure for gauge pressure
        virtual scalar pRefValue() const
        {
            return pRef_;
        }

        //- Absolute pressure (p + pRef) [Pa]
        virtual tmp<volScalarField> pAbs() const
        {
            if (pRefFound_)
            {
                return p() + pRef();
            }
            else
            {
                // Wrap a reference to original p field for efficiency
                return tmp<volScalarField>(p());
            }
        }

        //- Enthalpy/Internal energy/Temperature field
        //  (whichever is the input field to thermo)
        virtual const volScalarField& heT() const
        {
            return calculatesTFromhe() ? he() : T();
        }

        //- Return reference temperature to be added to T field
        virtual dimensionedScalar TRef() const
        {
            return dimensionedScalar("TRef", T().dimensions(), TRef_);
        }

        virtual scalar TRefValue() const
        {
            return TRef_;
        }

        //- Temperature
        virtual const volScalarField& T() const;

        //- Absolute temperature (T + TRef) [K]
        virtual tmp<volScalarField> TAbs() const
        {
            if (TRefFound_ && TRefValue() != scalar(0))
            {
                return T() + TRef();
            }
            else
            {
                // Wrap a reference to original T field for efficiency
                return tmp<volScalarField>(T());
            }
        }

        //- Temperature relative to TRef
        //  Non-const access allowed for transport equations
        virtual volScalarField& T();

        //- Properties dictionary
        virtual const IOdictionary& properties() const
        {
            return *this;
        }

        //- Access the dictionary
        virtual IOdictionary& properties()
        {
            return *this;
        }

        //- Thermal conductivity of mixture [W/m/K]
        virtual const volScalarField& kappa() const;


        //- Read thermophysical properties dictionary
        virtual bool read();


    // Member Operators

        //- Return material tables
        virtual materialTables& materials() const
        {
            materialPropertiesRequired();
            NotImplemented;
        }

        //- Return calculated model
        virtual tmp<volScalarField> operator[](const word& modelName) const;

        //- Disallow default bitwise assignment
        void operator=(const implementation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "basicThermo/basicThermoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
