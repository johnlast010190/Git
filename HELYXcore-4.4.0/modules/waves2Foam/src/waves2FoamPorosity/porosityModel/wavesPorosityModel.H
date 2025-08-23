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
    (c) held by original author
    (c) 2024 Engys Ltd.

Class
    Foam::wavesPorosityModel

Description


SourceFiles
    wavesPorosityModel.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef wavesPorosityModel_H
#define wavesPorosityModel_H

#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "cfdTools/general/include/fvCFD.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class wavesPorosityModel Declaration
\*---------------------------------------------------------------------------*/

class wavesPorosityModel
{
protected:

    // Protected data

        mutable volScalarField porosity_;

public:

    //- Runtime type information
    TypeName("wavesPorosityModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            wavesPorosityModel,
            wavesPorosityModel,
            (
                const fvMesh& mesh
            ),
            (mesh)
        );


    // Constructors

        //- Construct from components
        wavesPorosityModel
        (
            const fvMesh&
        );

        //- Disallow default bitwise copy construction
        wavesPorosityModel(const wavesPorosityModel&) = delete;


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<wavesPorosityModel> New
        (
            const fvMesh&
        );


    //- Destructor
    virtual ~wavesPorosityModel();


    // Member Functions

        //- The templated part had to be skipped, because it should also be a
        //  virtual function. I could not get this part working.
        //  NGJ, 29.10.2014

        //- mirror fvm::ddt with porosity
        virtual tmp<fvMatrix<vector>> ddt
        (
            VolField<vector>&
        ) = 0;

        //- mirror fvm::ddt with porosity
        virtual tmp<fvMatrix<vector>> ddt
        (
            const geometricOneField&,
            VolField<vector>&
        ) = 0;

        //- mirror fvm::ddt with porosity
        virtual tmp<fvMatrix<vector>> ddt
        (
            const dimensionedScalar&,
            VolField<vector>&
        ) = 0;

        //- mirror fvm::ddt with porosity
        virtual tmp<fvMatrix<vector>> ddt
        (
            const volScalarField&,
            VolField<vector>&
        ) = 0;

        //- Update the porosity field
        virtual void updatePorosity() = 0;

        //- Return a volScalarField with the porosities
        virtual const volScalarField& porosity() const = 0;

        //- Add the viscous and inertial resistance force contribution
        //  to the momentum equation
        virtual void addResistance(fvVectorMatrix& UEqn) const = 0;

        //- Add the viscous and inertial resistance force contribution
        //  to the tensorial diagonal
        virtual void addResistance
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU
        ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const wavesPorosityModel&) = delete;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
