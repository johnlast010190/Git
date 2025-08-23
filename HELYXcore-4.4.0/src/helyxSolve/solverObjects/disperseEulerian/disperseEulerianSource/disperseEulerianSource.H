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
    (c) 2015-2017 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::fv::disperseEulerianSource

Group
    grpFvOptionsSources

Description
    Calculates and applies the heat/mass/momentum transfer from disperse
    particles into the energy/mass/momentum equation of the continuous phase.

Usage
    Example usage:
    \verbatim
    fields          (T U alpha.<name>);   // List of field names
    \endverbatim

SourceFiles
    disperseEulerianSource.C

\*---------------------------------------------------------------------------*/

#ifndef disperseEulerianSource_H
#define disperseEulerianSource_H

#include "cfdTools/general/fvOptions/fvOption.H"
#include "solverObjects/disperseEulerian/disperseEulerian.H"
#include "solverObjects/disperseEulerian/phase/phase.H"
#include "solverObjects/disperseEulerian/disperseEulerianSource/heatTransfer/heatTransferModel/heatTransferModel.H"
#include "solverObjects/disperseEulerian/disperseEulerianSource/massTransfer/massTransferModel/massTransferModel.H"

#include "solverObjects/disperseEulerian/forceModels/dragModels/dragModel/dragModel.H"
#include "solverObjects/disperseEulerian/forceModels/liftModels/liftModel/liftModel.H"
#include "solverObjects/disperseEulerian/forceModels/turbulentDispersionModels/turbulentDispersionModel/turbulentDispersionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class disperseEulerianSource Declaration
\*---------------------------------------------------------------------------*/

class disperseEulerianSource
:
    public option
{
    // Private data

        //- Access to disperseEulerian
        const disperseEulerian& des_;

        //- Heat transfer model
        PtrList<decoupledEulerian::heatTransferModel> heatTrans_;

        //- mass transfer model
        PtrList<decoupledEulerian::massTransferModel> massTrans_;

        //- drag force model
        const PtrList<decoupledEulerian::dragModel> &forcesD_;

        //- lift force model
        const PtrList<decoupledEulerian::liftModel>& forcesL_;

        //- turbulent dispersion force model
        const PtrList<decoupledEulerian::turbulentDispersionModel>& forcesTD_;

        //- consider forces in boundary cells
        bool disableBoundaryForces_;

        //- list of patches at which to disable forces
        wordList patches_;


public:

    //- Runtime type information
    TypeName("disperseEulerianSource");


    // Constructors

        //- Construct from explicit source name and mesh
        disperseEulerianSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        disperseEulerianSource(const disperseEulerianSource&) = delete;


    // Member Functions

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames);

        //- Add explicit contribution to mass and temperature equations
        virtual void addSup
        (
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to momentum equation
        virtual void addSup
        (
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to compressible momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const disperseEulerianSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
