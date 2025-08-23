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
    (c) 2013-2019 OpenFOAM Foundation
    (c) 2022-2025 Engys Ltd.

Class
    Foam::pyrolysisChemistryModel

Description
    Pyrolysis chemistry model. It includes gas phase in the solid
    reaction.

SourceFiles
    pyrolysisChemistryModelI.H
    pyrolysisChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef pyrolysisChemistryModel_H
#define pyrolysisChemistryModel_H

#include "fields/volFields/volFieldsFwd.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"
#include "solidChemistryModel/solidChemistryModel.H"
#include "basicSolidChemistryModel/basicSolidChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                   Class pyrolysisChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class SolidThermo, class GasThermo>
class pyrolysisChemistryModel
:
    public solidChemistryModel<SolidThermo>
{
protected:

        //- List of gas species present in reaction system
        speciesTable pyrolisisGases_;

        //- Number of gas species
        label nGases_;

        //- Number of components being solved by ODE
        label nSpecie_;

        //- List of reaction rate per gas [kg/m3/s]
        PtrList<volScalarField::Internal> RRg_;


    // Protected Member Functions

        //- Write access to source terms for gases
        inline PtrList<volScalarField::Internal>& RRg();


private:

        //- List of accumulative solid concentrations
        mutable PtrList<volScalarField> Ys0_;


public:

    //- Runtime type information
    TypeName("pyrolysis");


    // Constructors

        //- Construct from thermo
        pyrolysisChemistryModel(const solidMulticomponentThermo& thermo);


    //- Destructor
    virtual ~pyrolysisChemistryModel();


    // Member Functions

        //- The number of solids
        inline label nSpecie() const;

        //- The number of solids
        inline label nGases() const;


        //- dc/dt = omega, rate of change in concentration, for each species
        virtual scalarField omega
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li,
            const bool updateC0 = false
        ) const;

        //- Return the reaction rate for reaction r
        virtual scalar omega
        (
            const Reaction<SolidThermo>& r,
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li,
            scalar& omegaf,
            scalar& omegar
        ) const;


        //- Calculates the reaction rates
        virtual void calculate();


        // Chemistry model functions

            //- Return const access to the chemical source terms for gases
            inline const volScalarField::Internal& RRg(const label i) const;

            //- Return total gas source term
            inline tmp<volScalarField::Internal> RRg() const;

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT);


        // ODE functions (overriding abstract functions in ODE.H)

            //- Number of ODE's to solve
            virtual label nEqns() const;

            virtual void derivatives
            (
                const scalar t,
                const scalarField& c,
                const label li,
                scalarField& dcdt
            ) const;

            virtual void jacobian
            (
                const scalar t,
                const scalarField& c,
                const label li,
                scalarField& dcdt,
                scalarSquareMatrix& dfdc
            ) const;

            virtual void solve
            (
                scalar& p,
                scalar& T,
                scalarField& c,
                const label li,
                scalar& deltaT,
                scalar& subDeltaT
            ) const
            {
                NotImplemented;
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pyrolysisChemistryModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pyrolysisChemistryModel/pyrolysisChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pyrolysisChemistryModel/pyrolysisChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
