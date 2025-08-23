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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2024 Engys Ltd

Class
    Foam::XiModels::algebraic

Description
    Simple algebraic model for Xi based on Gulders correlation
    with a linear correction function to give a plausible profile for Xi.
    See report TR/HGW/10 for details on the Weller two equations model.
    See \link XiModel.H \endlink for more details on flame wrinkling modelling.

SourceFiles
    algebraic.C

\*---------------------------------------------------------------------------*/

#ifndef algebraic_H
#define algebraic_H

#include "XiModel/XiModel.H"
#include "XiEqModel/XiEqModel.H"
#include "XiGModel/XiGModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiModels
{

/*---------------------------------------------------------------------------*\
                          Class algebraic Declaration
\*---------------------------------------------------------------------------*/

class algebraic
:
    public XiModel
{
    // Private data

        scalar XiShapeCoeff_;

        autoPtr<XiEqModel> XiEqModel_;
        autoPtr<XiGModel> XiGModel_;


public:

    //- Runtime type information
    TypeName("algebraic");


    // Constructors

        //- Construct from components
        algebraic
        (
            const dictionary& XiProperties,
            const psiuMulticomponentThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su,
            const volScalarField& rho,
            const volScalarField& b,
            const surfaceScalarField& phi
        );

        //- Disallow default bitwise copy construction
        algebraic(const algebraic&) = delete;


    //- Destructor
    virtual ~algebraic();


    // Member Functions

        //- Return the flame diffusivity
        virtual tmp<volScalarField> Db() const;

        //- Correct the flame-wrinking Xi
        virtual void correct();

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiProperties);

        //- Write fields of the XiEq model
        virtual void writeFields()
        {
            XiEqModel_().writeFields();
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const algebraic&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
