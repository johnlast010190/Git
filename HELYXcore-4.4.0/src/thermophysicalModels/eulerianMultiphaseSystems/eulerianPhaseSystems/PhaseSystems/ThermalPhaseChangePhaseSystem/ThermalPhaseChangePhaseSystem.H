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
    (c) 2015-2016 OpenFOAM Foundation

Class
    Foam::ThermalPhaseChangePhaseSystem

Description
    Class to provide interfacial heat and mass transfer between a number of
    phases according the interfacial temperature approximated by the saturation
    temperature.

    Currently only a single specified specie is considered volatile and changes
    phase, all other species are considered nonvolatile and do not
    affect the mass-transfer.

SourceFiles
    ThermalPhaseChangePhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef ThermalPhaseChangePhaseSystem_H
#define ThermalPhaseChangePhaseSystem_H

#include "../HeatAndMassTransferPhaseSystem/HeatAndMassTransferPhaseSystem.H"
#include "../../interfacialCompositionModels/eulerianSaturationModels/eulerianSaturationModel/eulerianSaturationModel.H"
#include "primitives/bools/Switch/Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class ThermalPhaseChangePhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class ThermalPhaseChangePhaseSystem
:
    public HeatAndMassTransferPhaseSystem<BasePhaseSystem>
{

protected:

    // Protected data

        //- Name of the volatile specie
        word volatile_;

        //- The saturation model used to evaluate Tsat = Tf
        autoPtr<eulerianSaturationModel> saturationModel_;

        // Mass transfer enabled
        Switch massTransfer_;

        //- Interfacial Mass transfer rate
        HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>
            iDmdt_;


public:

    // Constructors

        //- Construct from fvMesh
        ThermalPhaseChangePhaseSystem
        (
            const fvMesh&,
            const dictionary&
        );


    //- Destructor
    virtual ~ThermalPhaseChangePhaseSystem();


    // Member Functions

        //- Return the eulerianSaturationModel
        const eulerianSaturationModel& saturation() const;

        //- Return the heat transfer matrices
        virtual autoPtr<eulerianPhaseSystem::heatTransferTable> heatTransfer() const;

        //- Return the mass transfer matrices
        virtual autoPtr<eulerianPhaseSystem::massTransferTable> massTransfer() const;

        //- Return the interfacial mass flow rate
        virtual tmp<volScalarField> iDmdt(const phasePairKey& key) const;

        //- Return the total interfacial mass transfer rate for phase
        virtual tmp<volScalarField> iDmdt(const eulerianPhaseModel& phase) const;

        //- Correct the thermodynamics
        virtual void correctThermo();

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ThermalPhaseChangePhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
