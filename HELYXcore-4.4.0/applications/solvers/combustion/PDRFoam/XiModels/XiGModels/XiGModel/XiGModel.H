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
    (c) 2024 Engys Ltd

Class
    Foam::XiGModel

Description
    Base-class for all Xi generation models used by the b-Xi combustion model.
    See Technical Report SH/RE/01R for details on the PDR modelling. For details
    on the use of XiGModel see \link XiModel.H \endlink. The model available is
    \link instabilityG.H \endlink

SourceFiles
    XiGModel.C

\*---------------------------------------------------------------------------*/

#ifndef XiGModel_H
#define XiGModel_H

#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "psiuMulticomponentThermo/psiuMulticomponentThermo.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class XiGModel Declaration
\*---------------------------------------------------------------------------*/

class XiGModel
{

protected:

    // Protected data

        dictionary XiGModelCoeffs_;

        const psiuMulticomponentThermo& thermo_;
        const compressible::RASModel& turbulence_;
        const volScalarField& Su_;


public:

    //- Runtime type information
    TypeName("XiGModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            XiGModel,
            dictionary,
            (
                const dictionary& XiGProperties,
                const psiuMulticomponentThermo& thermo,
                const compressible::RASModel& turbulence,
                const volScalarField& Su
            ),
            (
                XiGProperties,
                thermo,
                turbulence,
                Su
            )
        );


    // Selectors

        //- Return a reference to the selected XiG model
        static autoPtr<XiGModel> New
        (
            const dictionary& XiGProperties,
            const psiuMulticomponentThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );


    // Constructors

        //- Construct from components
        XiGModel
        (
            const dictionary& XiGProperties,
            const psiuMulticomponentThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );

        //- Disallow default bitwise copy construction
        XiGModel(const XiGModel&) = delete;


    //- Destructor
    virtual ~XiGModel();


    // Member Functions

        //- Return the flame-wrinking genration rate
        virtual tmp<volScalarField> G() const = 0;

        //- Return the flame diffusivity
        virtual tmp<volScalarField> Db() const
        {
            return turbulence_.muEff();
        }

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiGProperties) = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const XiGModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
