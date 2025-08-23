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
    (c) 2019-2024 Engys Ltd.

Class
    Foam::fv::jouleHeatingSource

Group
    grpFvOptionsSources

Description
    Provides a Joule heating source according to:

    Differential form of Joule heating - power per unit volume:

    \f[
        \frac{d(P)}{d(V)} = J \cdot E
    \f]

    where \f$ J \f$ is the current density and \f$ E \f$ the electric field.
    If no magnetic field is present:

    \f[
        J = \sigma E
    \f]

    The electric field given by

    \f[
        E = \grad V
    \f]

    Therefore:

    \f[
        \frac{d(P)}{d(V)} = J \cdot E
                          = (sigma E) \cdot E
                          = (sigma \grad V) \cdot \grad V
    \f]


Usage
    (An)isotropic electrical conductivity and electrical potential are looked up
    - see electricalPotentialSolver
    \verbatim
    jouleHeatingSourceCoeffs
    {
    }
    \endverbatim

SourceFiles
    jouleHeatingSource.C

SeeAlso
    Foam::fv::electrostaticPotentialSolver

\*---------------------------------------------------------------------------*/

#ifndef fv_jouleHeatingSource_H
#define fv_jouleHeatingSource_H

#include "solverObjects/solverObject/solverObject.H"
#include "cfdTools/general/fvOptions/fvOption.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class jouleHeatingSource Declaration
\*---------------------------------------------------------------------------*/

class jouleHeatingSource
:
    public solverObject
{
public:

    //- Runtime type information
    TypeName("jouleHeatingSource");


    // Constructors

        //- Construct from explicit source name and mesh
        jouleHeatingSource
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        jouleHeatingSource(const jouleHeatingSource&) = delete;


    //- Destructor
    virtual ~jouleHeatingSource();


    // Member Functions

        //- Initialise after all fvOptions created
        virtual bool initialise() override;


        // Evaluate

            virtual void getSourceGraph
            (
                wordList& fieldNames,
                HashTable<wordList>& sourceDependencies
            ) override;

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            ) override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const jouleHeatingSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
