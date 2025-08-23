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
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::NoHeatTransfer

Group
    grpLagrangianIntermediateHeatTransferSubModels

Description
    Dummy heat transfer model for 'none'

\*---------------------------------------------------------------------------*/

#ifndef NoHeatTransfer_H
#define NoHeatTransfer_H

#include "submodels/Thermodynamic/HeatTransferModel/HeatTransferModel/HeatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class NoHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoHeatTransfer
:
    public HeatTransferModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from dictionary
        NoHeatTransfer(const dictionary&, CloudType& owner);

        //- Construct copy
        NoHeatTransfer(const NoHeatTransfer<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<HeatTransferModel<CloudType>> clone() const
        {
            return autoPtr<HeatTransferModel<CloudType>>
            (
                new NoHeatTransfer<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoHeatTransfer();


    // Member Functions

        //- Flag to indicate whether model activates heat transfer model
        virtual bool active() const;

        //- Nusselt number
        virtual scalar Nu(const scalar, const scalar) const;

        //- Prandtl number
        virtual scalar Pr() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "submodels/Thermodynamic/HeatTransferModel/NoHeatTransfer/NoHeatTransfer.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
