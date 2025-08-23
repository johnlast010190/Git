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
    (c) 2016-2017 OpenCFD Ltd.
    (c) 2015 IH-Cantabria

Class
    Foam::waveModels::regularWaveModel

Description

\*---------------------------------------------------------------------------*/

#ifndef waveModels_reguarWaveModel_H
#define waveModels_reguarWaveModel_H

#include "waveGenerationModels/base/waveGenerationModel/waveGenerationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{

/*---------------------------------------------------------------------------*\
                       Class regularWaveModel Declaration
\*---------------------------------------------------------------------------*/

class regularWaveModel
:
    public waveGenerationModel
{
protected:

    // Protected data

        //- Ramp time
        scalar rampTime_;

        //- Wave period
        scalar wavePeriod_;

        //- Wavelength
        scalar waveLength_;

        //- Wave phase
        scalar wavePhase_;


    // Protected Member Functions

        //- Return the time scaling coefficient
        virtual scalar timeCoeff(const scalar t) const;

        //- Return word description of wave type
        virtual word waveType() const;


public:

    //- Runtime type information
    TypeName("regularWaveModel");


    // Constructor

        //- Construct from components
        regularWaveModel
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const polyPatch& patch,
            const bool readFields = true
        );

        //- Disallow default bitwise copy construction
        regularWaveModel(const regularWaveModel&) = delete;


    //- Destructor
    virtual ~regularWaveModel();


    // Public Member Functions

        //- Read from dictionary
        virtual bool readDict(const dictionary& overrideDict);

        //- Info
        virtual void info(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const regularWaveModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
