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
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::massTransferModels::RanzMarshallSh

Description
    Ranz-Marshall correlation for turbulent heat transfer from the surface of a
    sphere to the surrounding fluid.

SourceFiles
    RanzMarshallSh.C

\*---------------------------------------------------------------------------*/

#ifndef RanzMarshallSh_H
#define RanzMarshallSh_H

#include "solverObjects/disperseEulerian/disperseEulerianSource/massTransfer/massTransferModel/massTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decoupledEulerian
{

class phase;

namespace massTransferModels
{

/*---------------------------------------------------------------------------*\
                          Class RanzMarshallSh Declaration
\*---------------------------------------------------------------------------*/

class RanzMarshallSh
:
    public massTransferModel
{
public:

    //- Runtime type information
    TypeName("RanzMarshall");


    // Constructors

        //- Construct from components
        RanzMarshallSh
        (
            const dictionary& dict,
            const phase& phase
        );


    //- Destructor
    virtual ~RanzMarshallSh();


    // Member Functions

        //- Compute the Sherwood number
        virtual tmp<volScalarField> Sh() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace massTransferModels
} // End namespace decoupledEulerian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
