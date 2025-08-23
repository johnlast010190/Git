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
    (c) 2024 Engys Ltd.

Class
    Foam::fv::acousticSource

Group
    grpFvOptionsSources

Description
    Simple volumetric source for generating acoustic noise for Perturbed
    Convective Wave Equation (PCWE)

    Example usage:
    \verbatim
    acousticSource1
    {
        type            acousticSource;
        active          yes;

        acousticSourceCoeffs
        {
            selectionMode      cellZone;
            cellZone           porosity;
            volumeMode         absolute;  // specific
            acousticSourceMode sinusoidal; // sinusoidalSuperposition  // table
            fieldName          psiAco;
            amplitude          1000;
            frequency          1000;
        }
    }
    \endverbatim


Note

SourceFiles
    acousticSource.C

\*---------------------------------------------------------------------------*/

#ifndef acousticSource_H
#define acousticSource_H

#include "cellSetOption/cellSetOption.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
              Class acousticSource Declaration
\*---------------------------------------------------------------------------*/

class acousticSource
:
    public cellSetOption
{

public:

    // Public data

        //- Enumeration for volume types
        enum volumeModeType
        {
            vmAbsolute,
            vmSpecific
        };

        //- Word list of volume mode type names
        static const NamedEnum<volumeModeType, 2> volumeModeTypeNames_;

        //- Enumeration for acoustic source mode types
        enum acousticSourceModeType
        {
            sinusoidal,
            sinusoidalSuperposition,
            table
        };

        //- Word list of acoustic source mode type names
        static const NamedEnum<acousticSourceModeType, 3> acousticSourceModeTypeNames_;

protected:

    // Protected data

        //- Volume mode
        volumeModeType volumeMode_;

        //- Volume normalisation
        scalar VDash_;

        //- Acoustic source mode
        acousticSourceModeType acousticSourceMode_;

        //- Acoustic source
        autoPtr<Function1<scalar>> AF1_;
        scalar A_;

        //- Amplitude of the acoustic source
        scalar amplitude_;

        //- Frequency of the acoustic source [Hz]
        scalar frequency_;

        //- Duration of the acoustic source
        scalar duration_;

        //- Minimum frequency for sinusoidal superposition of acoustic source
        scalar fMin_;

        //- Maximum frequency for sinusoidal superposition of acoustic source
        scalar fMax_;

        //- Number of components for sinusoidal superposition of acoustic source
        label nComps_;


    // Protected functions

        //- Set the acoustic source
        void updateSource(const fvMatrix<scalar>& eqn);

public:

    //- Runtime type information
    TypeName("acousticSource");


    // Constructors

        //- Construct from components
        acousticSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~acousticSource()
    {}


    // Member Functions

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames);


        // Evaluation

            //- Add contribution to equation
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldI
            );

            //- Add contribution to compressible equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldI
            );


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
