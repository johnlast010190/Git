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
    (c) 2017 IH-Cantabria
    (c) 2017 OpenCFD Ltd.

Class
    Foam::waveModels::Grimshaw

Description
    Grimshaw wave model

\*---------------------------------------------------------------------------*/

#ifndef waveModels_Grimshaw_H
#define waveModels_Grimshaw_H

#include "waveGenerationModels/base/solitaryWaveModel/solitaryWaveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{

/*---------------------------------------------------------------------------*\
                         Class Grimshaw Declaration
\*---------------------------------------------------------------------------*/

class Grimshaw
:
    public solitaryWaveModel
{
protected:

    // Protected Member Functions

        //- Wave height
        virtual scalar eta
        (
            const scalar H,
            const scalar h,
            const scalar x,
            const scalar y,
            const scalar theta,
            const scalar t,
            const scalar X0
        ) const;

        //-
        virtual scalar alfa
        (
            const scalar H,
            const scalar h
        ) const;

        //- Wave velocity
        virtual vector Uf
        (
            const scalar H,
            const scalar h,
            const scalar x,
            const scalar y,
            const scalar theta,
            const scalar t,
            const scalar X0,
            const scalar z
        ) const;

        //- Set the water level
        virtual void setLevel
        (
            const scalar t,
            const scalar tCoeff,
            scalarField& level
        ) const;

        //- Calculate the wave model velocity
        virtual void setVelocity
        (
            const scalar t,
            const scalar tCoeff,
            const scalarField& level
        );


public:

    //- Runtime type information
    TypeName("Grimshaw");

    //- Constructor
    Grimshaw
    (
        const dictionary& dict,
        const fvMesh& mesh,
        const polyPatch& patch,
        const bool readFields = true
    );

    //- Destructor
    virtual ~Grimshaw();


    // Public Member Functions

        //- Read from dictionary
        virtual bool readDict(const dictionary& overrideDict);

        //- Info
        virtual void info(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
