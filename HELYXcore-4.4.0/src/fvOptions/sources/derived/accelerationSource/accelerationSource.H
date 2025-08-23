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
    (c) 2018-2019 OpenFOAM Foundation

Class
    Foam::fv::accelerationSource

Description
    This fvOption applies an explicit acceleration force to components of the
    velocity field.

Usage
    Example usage:
    \verbatim
    accelerationSource
    {
        type        accelerationSource;
        active      on;
        selectionMode all;
        U           U;
        velocity    scale;
        //or acceleration scale;
        value       (-2.572 0 0);
        scale
        {
            type        halfCosineRamp;
            start       0;
            duration    10;
        }
    }
    \endverbatim

SourceFiles
    accelerationSource.C

\*---------------------------------------------------------------------------*/

#ifndef accelerationSource_H
#define accelerationSource_H

#include "cellSetOption/cellSetOption.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class accelerationSource Declaration
\*---------------------------------------------------------------------------*/

class accelerationSource
:
    public cellSetOption
{
    // Private Data

        //- Time-varying velocity
        autoPtr<Function1<vector>> velOrAccel_;


    // Private Member Functions

        //- Source term to momentum equation
        template<class AlphaRhoFieldType>
        void add
        (
            const AlphaRhoFieldType& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Source term to momentum equation
        template<class AlphaRhoFieldType>
        void add
        (
            const AlphaRhoFieldType& rho,
            fvBlockMatrix<vector>& eqn,
            const label fieldi
        );


public:

    //- Runtime type information
    TypeName("accelerationSource");


    // Constructors

        //- Construct from components
        accelerationSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~accelerationSource()
    {}


    // Member Functions

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames);


         // Add explicit and implicit contributions

            //- Source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Source term to phase momentum equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Source term to momentum equation
            virtual void addSup
            (
                fvBlockMatrix<vector>& eqn,
                const label fieldi
            );

            //- Source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvBlockMatrix<vector>& eqn,
                const label fieldi
            );


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "accelerationSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
