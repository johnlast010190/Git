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
    (c) 2010-2024 Engys Ltd.
    (c) 1991-2008 OpenCFD Ltd.

Class
    Foam::concentrationTransport

Description
    A function object for solving passive scalar transport. Where applicable
    optional entries include their default values. Some optional entries have
    no or NULL defaults

    Example of the function object specification:
    @verbatim
    tracer
    {
        type            concentrationTransport;
        region          region0;

        fieldName       Trace;
        phiName         phi; //optional
        rhoName         rho; //optional
    }
    @endverbatim

SourceFiles
    concentrationTransport.C

\*---------------------------------------------------------------------------*/

#ifndef concentrationTransport_H
#define concentrationTransport_H

#include "solverObjects/solverObject/solverObject.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "cfdTools/general/include/fvCFD.H"
#include "sets/topoSets/cellSet.H"
#include "cfdTools/general/fvOptions/fvOptionList.H"
#include "dynamicFvMesh/dynamicFvMesh.H"
#include "primitives/strings/lists/hashedWordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class concentrationTransport Declaration
\*---------------------------------------------------------------------------*/

class concentrationTransport
:
    public solverObject
{
protected:

    // Protected data

        //- Name of the transport field
        word fieldName_;

        //- Name of flux field (optional)
        word phiName_;

        //- Name of density field (optional)
        word rhoName_;

        //- laminar diffusion
        autoPtr<dimensionedScalar> D_;

        //- turbulent Schmidt number
        autoPtr<dimensionedScalar> Sct_;

        //- Schmidt number - optional
        autoPtr<dimensionedScalar> Sc_;

        //- phase name - for multiphase
        autoPtr<word> phaseName_;

        //- phase fraction field - for multiphase
        const volScalarField* phaseFieldPtr_;

        //- Modify deltaT to speed-up passive transport (optional)
        label solveSpeedup_;

        //- Allow sub-cycling with specified number of time steps
        label nSubCycles_;

        //- Stabilize phasic transport (optional)
        dimensionedScalar residualPhaseField_;

        //- Specify list of boundaries with zero diffusion flux
        hashedWordList noDiffusionBoundaries_;

        //- use old-style phasic transport (more stable, less accurate) y/n
        bool stablePhasic_;

        //- use phasic fvOption sources y/n
        bool phasicSources_;


    // Private Member Functions

        const volScalarField& phase() const;


public:

    //- Runtime type information
    TypeName("concentrationTransport");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        concentrationTransport
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary&
        );


    //- Destructor
    virtual ~concentrationTransport();


    // Member Functions

         //- Returns diffusivity
        tmp<volScalarField> diffusivity() const;

        virtual bool initialise() override;

        //- Identify any solves handled by this solverObject, their levels of
        // iteration and dependent solves
        virtual void getSolveGraph
        (
            wordList& solveNames,
            HashTable<wordList>& requiredDependencies,
            HashTable<wordList>& optionalDependencies,
            HashTable<wordList>& correctorMembers
        ) override;

        //- Update solver
        virtual void correct(const word&, const word&) override;

        //- assemble and solve inisde separate function since
        //  optional sub-cycling otherwise duplicates all code
        virtual void solveEquation
        (
            volScalarField& field,
            const surfaceScalarField& phi,
            fv::options& fvOptions
        );

        //- Sample and write
        virtual void write() override;

        //- Read the parameters
        virtual void read(const dictionary&) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
