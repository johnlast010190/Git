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
    Foam::porosityZones

Description
    A simplified version of the original OF porosity module, now added the
    functionality of runTime selection of the type of resistance coefficients.

SourceFiles
    porosityZones.C

\*---------------------------------------------------------------------------*/

#ifndef porosityZones_H
#define porosityZones_H

#include "porosityZone.H"
#include "db/IOobjects/IOPtrList/IOPtrList.H"

#include "cfdTools/general/include/fvCFD.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class porosityZones Declaration
\*---------------------------------------------------------------------------*/

class porosityZones
:
    public IOPtrList<porosityZone>
{
    // Private data

        //- Reference to the finite volume mesh this zone is part of
        const fvMesh& mesh_;


    // Private Member Functions

        //- Modify time derivative elements
        template<class Type>
        void modifyDdt(fvMatrix<Type>&) const;


public:

    // Constructors

        //- Construct from fvMesh
        //  with automatically constructed coordinate systems list
        porosityZones(const fvMesh&);

        //- Disallow default bitwise copy construction
        porosityZones(const porosityZones&) = delete;


    // Member Functions

        //- mirror fvm::ddt with porosity
        template<class Type>
        tmp<fvMatrix<Type>> ddt
        (
            VolField<Type>&
        );

        //- mirror fvm::ddt with porosity
        template<class Type>
        tmp<fvMatrix<Type>> ddt
        (
            const geometricOneField&,
            VolField<Type>&
        );

        //- mirror fvm::ddt with porosity
        template<class Type>
        tmp<fvMatrix<Type>> ddt
        (
            const dimensionedScalar&,
            VolField<Type>&
        );

        //- mirror fvm::ddt with porosity
        template<class Type>
        tmp<fvMatrix<Type>> ddt
        (
            const volScalarField&,
            VolField<Type>&
        );

        //- Return a volScalarField with the porosities
        tmp<volScalarField> porosity() const;

        //- Add the viscous and inertial resistance force contribution
        //  to the momentum equation
        void addResistance(fvVectorMatrix& UEqn) const;

        //- Add the viscous and inertial resistance force contribution
        //  to the tensorial diagonal
        void addResistance
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU
        ) const;

        //- read modified data
        virtual bool readData(Istream&);

        //- write data
        bool writeData(Ostream&, bool subDict = true) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const porosityZones&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "porosityZonesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
