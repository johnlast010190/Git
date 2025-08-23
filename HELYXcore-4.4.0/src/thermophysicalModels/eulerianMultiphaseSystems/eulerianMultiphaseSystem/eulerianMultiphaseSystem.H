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
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd

Class
    Foam::eulerianMultiphaseSystem

Description
    Class which solves the volume fraction equations for two phases.

SourceFiles
    eulerianMultiphaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseSystem_H
#define multiphaseSystem_H

#include "../eulerianPhaseSystems/eulerianPhaseSystem/eulerianPhaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class eulerianDragModel;
class eulerianVirtualMassModel;

/*---------------------------------------------------------------------------*\
                  Class eulerianMultiphaseSystem Declaration
\*---------------------------------------------------------------------------*/

class eulerianMultiphaseSystem
:
    public eulerianPhaseSystem
{
    // Private member functions

        //- Return the drag coefficient for phase pair
        virtual tmp<volScalarField> Kd(const phasePairKey& key) const = 0;

        //- Return the face drag coefficient for phase pair
        virtual tmp<surfaceScalarField> Kdf(const phasePairKey& key) const = 0;

        //- Return the virtual mass coefficient for phase pair
        virtual tmp<volScalarField> Vm(const phasePairKey& key) const = 0;

        //- Return the face virtual mass coefficient for phase pair
        virtual tmp<surfaceScalarField> Vmf(const phasePairKey& key) const = 0;

        //- Return the turbulent diffusivity for phase pair
        //  Multiplies the phase-fraction gradient
        virtual tmp<volScalarField> D(const phasePairKey& key) const = 0;

        //- Return the interfacial mass flow rate for phase pair
        virtual tmp<volScalarField> dmdt(const phasePairKey& key) const = 0;


public:

    //- Runtime type information
    TypeName("multiphaseSystem");

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            eulerianMultiphaseSystem,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& coeffs
            ),
            (mesh, coeffs)
        );


    // Constructors

        //- Construct from fvMesh
        eulerianMultiphaseSystem(const fvMesh&, const dictionary&);


    //- Destructor
    virtual ~eulerianMultiphaseSystem();


    // Selectors

        static autoPtr<eulerianMultiphaseSystem> New
        (
            const fvMesh& mesh, const basicThermo& thermo
        );


    // Member Functions

        //- Return the drag coefficient for all phase-pairs
        virtual const eulerianPhaseSystem::KdTable& Kds() const = 0;

        //- Return the drag coefficient for phase
        virtual tmp<volScalarField> Kd(const eulerianPhaseModel& phase) const = 0;

        //- Return the combined force (lift + wall-lubrication) for phase pair
        virtual autoPtr<PtrList<Foam::volVectorField>> Fs() const = 0;

        //- Return the turbulent dispersion force on faces for phase pair
        virtual autoPtr<PtrList<Foam::surfaceScalarField>> phiDs
        (
            const PtrList<volScalarField>& rAUs
        ) const = 0;

        //- Return true if there is mass transfer for phase
        virtual bool transfersMass(const eulerianPhaseModel& phase) const = 0;

        //- Return the total interfacial mass transfer rate for phase
        virtual tmp<volScalarField> dmdt(const eulerianPhaseModel& phase) const = 0;

        //- Return the momentum transfer matrices
        virtual autoPtr<momentumTransferTable> momentumTransfer() const = 0;

        //- Return the heat transfer matrices
        virtual autoPtr<heatTransferTable> heatTransfer() const = 0;

        //- Return the mass transfer matrices
        virtual autoPtr<massTransferTable> massTransfer() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
