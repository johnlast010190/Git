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
    (c) 2011-2022 OpenFOAM Foundation

Class
    Foam::distributionMapLagrangian

Description
    Class containing mesh-to-mesh mapping information for particles

SourceFiles
    distributionMapLagrangian.C

\*---------------------------------------------------------------------------*/

#ifndef distributionMapLagrangian_H
#define distributionMapLagrangian_H

#include "meshes/polyMesh/polyDistributionMap/distributionMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class distributionMapLagrangian Declaration
\*---------------------------------------------------------------------------*/

class distributionMapLagrangian
{
    // Private data

        //- Map to distribute particles
        const distributionMap particleMap_;

        //- Per element in subsetted mesh the cell label
        const labelListList constructCellLabels_;


public:

    // Constructors

        //- Construct from components
        distributionMapLagrangian
        (
            const label nNewParticles,
            labelListList&& subParticleMap,
            labelListList&& constructParticleMap,
            labelListList&& constructCellLabels
        )
        :
            particleMap_(nNewParticles, subParticleMap, constructParticleMap),
            constructCellLabels_(constructCellLabels)
        {}


    // Member Functions

        // Access

            //- Distribution map
            const distributionMap& particleMap() const
            {
                return particleMap_;
            }

            //- Per received particle the destination cell label
            const labelListList& constructCellLabels() const
            {
                return constructCellLabels_;
            }


        // Edit

            //- Distribute list of lagrangian data
            template<class T>
            void distributeLagrangianData(List<T>& lst) const
            {
                particleMap_.distribute(lst);
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
