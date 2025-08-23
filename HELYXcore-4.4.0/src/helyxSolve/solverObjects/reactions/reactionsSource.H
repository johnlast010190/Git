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
    (c) 2013-2022 OpenFOAM Foundation
    (c) 2019-2024 Engys Ltd.

Class
    Foam::fv::reactionsSource

Group
    grpFvOptionsSources

Description
    Thermal and species source terms from combustion model.

SourceFiles
    reactionsSource.C

\*---------------------------------------------------------------------------*/

#ifndef reactionsSource_H
#define reactionsSource_H

#include "solverObjects/solverObject/solverObject.H"
#include "combustionModel/combustionModel.H"
#include "rhoMulticomponentThermo/rhoMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class reactionsSource Declaration
\*---------------------------------------------------------------------------*/

class reactionsSource
:
    public solverObject
{
protected:


    // Protected Member Data

        //- Reaction model
        autoPtr<combustionModel> reaction_;

        //- Pointer to solution control
        solutionControl* solnControlPtr_;

        //- Is it reactingFoam eqvivalent?
        Switch isReactingFoam_;


public:

    //- Runtime type information
    TypeName("reactionsSource");


    // Constructors

        //- Construct from components
        reactionsSource
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construct
        reactionsSource(reactionsSource&) = delete;


    //- Destructor
    virtual ~reactionsSource()
    {}


    // Solver Functions

        virtual bool initialise() override;

        virtual void read(const dictionary& dict) override;

        virtual void getSolveGraph
        (
            wordList& solveNames,
            HashTable<wordList>& requiredDependencies,
            HashTable<wordList>& optionalDependencies,
            HashTable<wordList>& correctorMembers
        ) override;

        virtual void correct(const word&, const word&) override;


    // Source term functions

        virtual void getSourceGraph
        (
            wordList& fields,
            HashTable<wordList>& sourceDependencies
        ) override;

        //- Add contribution to compressible equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldI
        ) override;

        //- Calculates the local deltaT for the flow solver
        virtual bool setRDeltaT() override;

        //- Is it reactingFoam eqvivalent?
        bool reactingFoam() const
        {
            return isReactingFoam_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(reactionsSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
