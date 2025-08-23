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
    (c) 2010-2024 Engys Ltd.
    (c) 2012-2016 OpenFOAM Foundation

Class
    Foam::fv::explicitPorositySource

Group
    grpFvOptionsSources

Description
    Explicit porosity source

Usage
    Example usage, here employing the Darcy-Forchheimer model:
    \verbatim
    explicitPorositySourceCoeffs
    {
        type            DarcyForchheimer;

        DarcyForchheimerCoeffs
        {
            d   d [0 -2 0 0 0 0 0] (5e7 -1000 -1000);
            f   f [0 -1 0 0 0 0 0] (0 0 0);
            referenceFrame  porousFrame;
        }
    }
    \endverbatim

Note:
    The porous region must be selected as a cellZone.

SourceFiles
    explicitPorositySource.C

\*---------------------------------------------------------------------------*/

#ifndef explicitPorositySource_H
#define explicitPorositySource_H

#include "cellSetOption/cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class porosityModel;

namespace fv
{


/*---------------------------------------------------------------------------*\
                   Class explicitPorositySource Declaration
\*---------------------------------------------------------------------------*/

class explicitPorositySource
:
    public cellSetOption
{
protected:

    // Protected data

        //- Run-time selectable porosity model
        autoPtr<porosityModel> porosityPtr_;

private:

        //- Adjoint velocity field prefix
        word adjointVelocityPrefix_;

        //- Velocity field name
        word UName_;


public:

    //- Runtime type information
    TypeName("explicitPorositySource");


    // Constructors

        //- Construct from components
        explicitPorositySource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        explicitPorositySource(const explicitPorositySource&) = delete;


    //- Destructor
    virtual ~explicitPorositySource()
    {}


    // Member Functions

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames);


        // Add explicit and implicit contributions

            //- Add implicit contribution to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Add implicit contribution to block Momentum equation
            virtual void addSup
            (
                fvBlockMatrix<vector>& eqn,
                const label fieldi
            );

            //- Add implicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Add implicit contribution to block Momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvBlockMatrix<vector>& eqn,
                const label fieldi
            );

            //- Add implicit contribution to phase momentum equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- always return true for adjoint velocity fields
            virtual label applyToField
            (
                const word& fieldName, const word& regionName
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const explicitPorositySource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
