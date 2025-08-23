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
    (c) 2012-2014 OpenFOAM Foundation

Class
    Foam::porosityModelList

Description
    List container for porosity models

SourceFiles
    porosityModelList.C

\*---------------------------------------------------------------------------*/

#ifndef porosityModelList_H
#define porosityModelList_H

#include "fvMesh/fvMesh.H"
#include "db/dictionary/dictionary.H"
#include "fvMatrices/fvMatricesFwd.H"
#include "fvMatrices/fvBlockMatrix/fvBlockMatrix.H"
#include "cfdTools/general/porosityModel/porosityModel/porosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class porosityModelList;
Ostream& operator<<(Ostream& os, const porosityModelList& models);

/*---------------------------------------------------------------------------*\
                      Class porosityModelList Declaration
\*---------------------------------------------------------------------------*/

class porosityModelList
:
    public PtrList<porosityModel>
{
protected:

    // Protected data
        const objectRegistry& obr_;

        //- Reference to the mesh database
        const fvMesh& mesh_;


public:


    // Constructors

        //- Constructor
        porosityModelList
        (
            const objectRegistry& obr,
            const fvMesh& mesh,
            const dictionary& dict,
            bool readFromFvOptions = false
        );

        //- Disallow default bitwise copy construction
        porosityModelList(const porosityModelList&) = delete;

    //- Destructor
    ~porosityModelList();


    // Member Functions

        //- Return active status
        bool active(const bool active = false) const;

        //- Reset the source list
        void reset(const dictionary& dict, bool readFromFvOptions);

        //- Add resistance
        void addResistance(fvVectorMatrix& UEqn);

        //- Add resistance
        void addResistance
        (
            fvBlockMatrix<vector>&
        );

        //- Add resistance
        void addResistance
        (
            fvVectorMatrix& UEqn,
            const volScalarField& rho,
            const volScalarField& mu
        );

        //- Add resistance
        void addResistance
        (
            fvBlockMatrix<vector>&,
            const volScalarField&,
            const volScalarField&
        );

        //- Add resistance
        void addResistance
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU,
            bool correctAUprocBC = true
        );

        //- Add adjoint resistance
        void addAdjointResistance
        (
            fvVectorMatrix& UaEqn,
            const volVectorField& Uprimal
        );

        //- Add adjoint resistance
        void addAdjointResistance
        (
            fvVectorMatrix& UaEqn,
            const volScalarField& rho,
            const volScalarField& mu,
            const volVectorField& Uprimal
        );

        //- Add adjoint resistance
        void addAdjointResistance
        (
            fvBlockMatrix<vector>&,
            const volVectorField&
        );


        //- Add adjoint resistance
        void addAdjointResistance
        (
            fvBlockMatrix<vector>& UaEqn,
            const volScalarField& rho,
            const volScalarField& mu,
            const volVectorField& Uprimal
        );

        //- Add forces and moments
        void addForceMoment
        (
            const volVectorField& U,
            const volScalarField& rho,
            const volScalarField& mu,
            vectorField& porousForce,
            vectorField& porousMoment
         );

        //- Find porous index given a name
        label findPorosityModelID
        (
            const word& name
        ) const;

        //- Return a list of porous models
        wordList names() const;

        //- return true if selection mode is cellSet
        bool selectionMode(const dictionary& dict);

        // I-O

            //- Read dictionary
            bool read(const dictionary& dict);

            //- Write data to Ostream
            bool writeData(Ostream& os) const;


    // Member Operators

        //- Return const and non-const reference to porosityModel by index.
        using PtrList<porosityModel>::operator[];

        //- Return const reference to ZoneType by name.
        const porosityModel& operator[](const word&) const;

        //- Ostream operator
        friend Ostream& operator<<
        (
            Ostream& os,
            const porosityModelList& models
        );

        //- Disallow default bitwise assignment
        void operator=(const porosityModelList&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
