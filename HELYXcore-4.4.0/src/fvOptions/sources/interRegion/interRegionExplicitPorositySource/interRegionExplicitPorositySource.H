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
    (c) 2012-2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::fv::interRegionExplicitPorositySource

Group
    grpFvOptionsSources

Description
    Inter-region explicit porosity source.

    Sources described by, for example using the DarcyForchheimer model:

    \verbatim
    interRegionExplicitPorositySourceCoeffs
    {
        type            DarcyForchheimer;
        DarcyForchheimerCoeffs
        {
            d   d [0 -2 0 0 0 0 0] (5e7 -1000 -1000);
            f   f [0 -1 0 0 0 0 0] (0 0 0);

            coordinateSystem
            {
                e1  (0.70710678 0.70710678 0);
                e2  (0 0 1);
            }
        }
    }
    \endverbatim

Note
    The porous region must be selected as a cellZone.

SourceFiles
    interRegionExplicitPorositySource.C

\*---------------------------------------------------------------------------*/

#ifndef interRegionExplicitPorositySource_H
#define interRegionExplicitPorositySource_H

#include "interRegionOption/interRegionOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class porosityModel;

namespace fv
{


/*---------------------------------------------------------------------------*\
                   Class interRegionExplicitPorositySource Declaration
\*---------------------------------------------------------------------------*/

class interRegionExplicitPorositySource
:
    public interRegionOption
{
protected:

    // Protected data

        //- Run-time selectable porosity model
        autoPtr<porosityModel> porosityPtr_;

        //- First iteration
        bool firstIter_;

        //- Velocity field name, default = U
        word UName_;

        //- Dynamic viscosity field name (compressible case only)
        //  default = mu
        word muName_;


    // Protected Member Functions

        //- Initialise
        void init();


public:

    //- Runtime type information
    TypeName("interRegionExplicitPorositySource");


    // Constructors

        //- Construct from components
        interRegionExplicitPorositySource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        interRegionExplicitPorositySource
        (
            const interRegionExplicitPorositySource&
        ) = delete;


    //- Destructor
    virtual ~interRegionExplicitPorositySource()
    {}


    // Member Functions

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames);

        //- Add explicit and implicit contributions
        virtual void addSup
        (
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Add explicit and implicit contributions to compressible equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Read dictionary
        virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const interRegionExplicitPorositySource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
