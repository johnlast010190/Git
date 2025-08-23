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
    (c) 2015-2017 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::fv::tabulatedAccelerationSource

Group
    grpFvOptionsSources

Description
    Solid-body 6-DoF acceleration source

Usage
    Example usage:
    \verbatim
    SBM
    {
        type            tabulatedAccelerationSource;
        active          yes;

        timeDataFileName "constant/acceleration-terms.dat";
    }
    \endverbatim

SourceFiles
    tabulatedAccelerationSource.C

\*---------------------------------------------------------------------------*/

#ifndef tabulatedAccelerationSource_H
#define tabulatedAccelerationSource_H

#include "cfdTools/general/fvOptions/fvOption.H"
#include "sources/derived/tabulatedAccelerationSource/tabulated6DoFAcceleration/tabulated6DoFAcceleration.H"
#include "dimensionedTypes/dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class tabulatedAccelerationSource Declaration
\*---------------------------------------------------------------------------*/

class tabulatedAccelerationSource
:
    public option
{
protected:

    // Protected data

        //- Run-time selectable acceleration model
        tabulated6DoFAcceleration motion_;

        //- Velocity field name, default = U
        word UName_;

        dimensionedVector g0_;


private:

    // Private Member Functions

        //- Source term to momentum equation
        template<class RhoFieldType>
        void addSup
        (
            const RhoFieldType& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );


public:

    //- Runtime type information
    TypeName("tabulatedAccelerationSource");


    // Constructors

        //- Construct from components
        tabulatedAccelerationSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        tabulatedAccelerationSource
        (
            const tabulatedAccelerationSource&
        ) = delete;


    //- Destructor
    virtual ~tabulatedAccelerationSource()
    {}


    // Member Functions

        //- Initialise after all fvOptions created
        virtual bool initialise() override;

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames) override;

        //- Source term to momentum equation
        virtual void addSup
        (
            fvMatrix<vector>& eqn,
            const label fieldi
        ) override;

        //- Source term to compressible momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        ) override;

        //- Read dictionary
        virtual bool read(const dictionary& dict) override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const tabulatedAccelerationSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "tabulatedAccelerationSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
