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
    (c) 2021-2025 Engys Ltd.

Class
    Foam::equationOfState

Description
    Class for selecting material models.

SourceFiles
    equationOfState.C

\*---------------------------------------------------------------------------*/

#ifndef equationOfState_H
#define equationOfState_H

#include "materialModels/materialTables/materialTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class equationOfState Declaration
\*---------------------------------------------------------------------------*/

class equationOfState
{

public:

    //- Runtime type information
    TypeName("equationOfState");

    //- Lists models associated with equation of state
    static const wordList listModels()
    {
        return
            wordList
            ({
                rhoModel::typeName,
                spContributionModel::typeName,
                hContributionModel::typeName,
                CpContributionModel::typeName,
                svContributionModel::typeName,
                eContributionModel::typeName,
                CvContributionModel::typeName,
                psiModel::typeName,
                ZModel::typeName,
                CpMCvModel::typeName,
                alphavModel::typeName
            });
    }

    //- Destructor
    virtual ~equationOfState();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
