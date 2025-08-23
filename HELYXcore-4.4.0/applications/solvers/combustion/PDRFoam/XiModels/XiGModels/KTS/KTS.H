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
    Foam::XiGModels::KTS

Description
    Simple Kolmogorov time-scale (KTS) model for the flame-wrinling generation
    rate.

SourceFiles
    KTS.C

\*---------------------------------------------------------------------------*/

#ifndef KTS_H
#define KTS_H

#include "XiGModel/XiGModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiGModels
{

/*---------------------------------------------------------------------------*\
                          Class KTS Declaration
\*---------------------------------------------------------------------------*/

class KTS
:
    public XiGModel
{
    // Private data

        scalar GEtaCoef_;


public:

    //- Runtime type information
    TypeName("KTS");


    // Constructors

        //- Construct from components
        KTS
        (
            const dictionary& XiGProperties,
            const psiuMulticomponentThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );

        //- Disallow default bitwise copy construction
        KTS(const KTS&) = delete;


    //- Destructor
    virtual ~KTS();


    // Member Functions

        //- Return the flame-wrinking generation rate
        virtual tmp<volScalarField> G() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiGProperties);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const KTS&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiGModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
