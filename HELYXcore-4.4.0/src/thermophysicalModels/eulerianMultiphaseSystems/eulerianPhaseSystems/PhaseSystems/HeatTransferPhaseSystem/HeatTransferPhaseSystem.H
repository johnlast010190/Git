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
    Foam::HeatTransferPhaseSystem

Description
    Class which models interfacial heat transfer between a number of phases.

SourceFiles
    HeatTransferPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef HeatTransferPhaseSystem_H
#define HeatTransferPhaseSystem_H

#include "../../eulerianPhaseSystem/eulerianPhaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class eulerianBlendingMethod;
template<class modelType> class EulerianBlendedInterfacialModel;
class eulerianHeatTransferModel;

/*---------------------------------------------------------------------------*\
                    Class HeatTransferPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class HeatTransferPhaseSystem
:
    public BasePhaseSystem
{
protected:

    // Protected typedefs

        typedef HashTable
        <
            autoPtr<EulerianBlendedInterfacialModel<eulerianHeatTransferModel>>,
            phasePairKey,
            phasePairKey::hash
        > heatTransferModelTable;


    // Protected data

        // Sub Models

            //- Heat transfer models
            heatTransferModelTable heatTransferModels_;


public:

    // Constructors

        //- Construct from fvMesh
        HeatTransferPhaseSystem(const fvMesh&, const dictionary&);


    //- Destructor
    virtual ~HeatTransferPhaseSystem();


    // Member Functions

        //- Return true if there is mass transfer for phase
        virtual bool transfersMass(const eulerianPhaseModel& phase) const;

        //- Return the interfacial mass flow rate
        virtual tmp<volScalarField> dmdt
        (
            const phasePairKey& key
        ) const;

        //- Return the total interfacial mass transfer rate for phase
        virtual tmp<volScalarField> dmdt(const eulerianPhaseModel& phase) const;

        //- Return the heat transfer matrices
        virtual autoPtr<eulerianPhaseSystem::heatTransferTable>
            heatTransfer() const;

        //- Return the mass transfer matrices
        virtual autoPtr<eulerianPhaseSystem::massTransferTable>
            massTransfer() const;

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HeatTransferPhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
