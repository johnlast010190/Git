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
    Foam::fv::turbulenceControlledDecaySource

Group
    grpFvOptionsSources

Description
    Calculates and applies controlled decay source term effectively
    setting a lower limit for turbulent decay by couter-acting
    dissipation with same term of opposite sign and ambient values.
    (Currently only supports k-omega turbulence models!)

Usage
    Example usage:
    \verbatim
    fields          (k omega);  // List of field names
    kName           k;          // name of tke field
    tName           omega;      // name of turb diss field
    kAmb [0 2 -2 0 0] 10;       // ambient value for k
    tAmb [0 0 -1 0 0] 1;        // ambient value for omega
    \endverbatim

SourceFiles
    turbulenceControlledDecaySource.C

\*---------------------------------------------------------------------------*/

#ifndef turbulenceControlledDecaySource_H
#define turbulenceControlledDecaySource_H

#include "cfdTools/general/fvOptions/fvOption.H"
#include "cellSetOption/cellSetOption.H"
#include "memory/autoPtr/autoPtr.H"
#include "sets/topoSetSource/topoSetSource.H"
#include "sets/topoSets/cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class turbulenceControlledDecaySource Declaration
\*---------------------------------------------------------------------------*/

class turbulenceControlledDecaySource
:
    public cellSetOption
{
    // Private data

        // k field name
        word kName_;

        // epsilon/omega field name
        word tName_;

        // k ambient
        dimensionedScalar kAmb_;

        // epsilon/omega ambient
        dimensionedScalar tAmb_;


public:

    //- Runtime type information
    TypeName("turbulenceControlledDecaySource");


    // Constructors

        //- Construct from explicit source name and mesh
        turbulenceControlledDecaySource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        turbulenceControlledDecaySource
        (
            const turbulenceControlledDecaySource&
        ) = delete;


    // Member Functions

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames);


        // Evaluate

            //- Add explicit contribution to turbulence equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            );
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldi
            );

        // I-O

            //- Read dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const turbulenceControlledDecaySource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
