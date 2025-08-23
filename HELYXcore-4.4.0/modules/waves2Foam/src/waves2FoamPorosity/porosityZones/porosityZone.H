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
    Foam::porosityZone

Description
    A simplified version of the original OF porosity module, now added the
    functionality of runTime selection of the type of resistance coefficients.

See Also
    porosityZones and coordinateSystems

SourceFiles
    porosityZone.C
    porosityZoneTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef porosityZone_H
#define porosityZone_H

#include "db/dictionary/dictionary.H"
#include "coordinate/systems/coordinateSystem.H"
#include "coordinate/systems/coordinateSystems.H"
#include "primitives/strings/lists/wordList.H"
#include "primitives/ints/lists/labelList.H"
#include "dimensionedTypes/dimensionedScalar/dimensionedScalar.H"
#include "dimensionedTypes/dimensionedTensor/dimensionedTensor.H"
#include "fields/Fields/primitiveFieldsFwd.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fvMatrices/fvMatricesFwd.H"

#include "fvMesh/fvMesh.H"

#include "porosityCoefficients/porosityCoefficient.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                        Class porosityZone Declaration
\*---------------------------------------------------------------------------*/

class porosityZone
{
    // Private data

        //- Name of this zone
        word name_;

        //- Reference to the finite volume mesh this zone is part of
        const fvMesh& mesh_;

        //- Dictionary containing the parameters
        dictionary dict_;

        //- Cell zone ID
        label cellZoneID_;

        //- Coordinate system used for the zone (Cartesian)
        coordinateSystem coordSys_;

        //- Porosity of the zone (0 < porosity <= 1)
        //  Placeholder for treatment of temporal terms.
        scalar porosity_;

        scalar addedMassCoeff_;

        //- Darcy coefficient
        dimensionedTensor D_;

        //- Forchheimer coefficient
        dimensionedTensor F_;


    // Private Member Functions

        //- adjust negative resistance values to be multiplier of max value
        static void checkNegativeResistance(dimensionedVector& resist);

        //- Viscous and inertial resistance
        template<class RhoFieldType>
        void addViscousInertialResistance
        (
            scalarField& Udiag,
            vectorField& Usource,
            const labelList& cells,
            const scalarField& V,
            const RhoFieldType& rho,
            const scalarField& mu,
            const vectorField& U
        ) const;

        //- Viscous and inertial resistance
        template<class RhoFieldType>
        void addViscousInertialResistance
        (
            tensorField& AU,
            const labelList& cells,
            const RhoFieldType& rho,
            const scalarField& mu,
            const vectorField& U
        ) const;


public:

    // Constructors

        //- Construct from components
        porosityZone(const word& name, const fvMesh&, const dictionary&);

        //- Disallow default bitwise copy construction
        porosityZone(const porosityZone&) = delete;

        //- Return clone
        autoPtr<porosityZone> clone() const
        {
            notImplemented("autoPtr<porosityZone> clone() const");
            return autoPtr<porosityZone>(nullptr);
        }

        //- Return pointer to new porosityZone created on freestore from Istream
        class iNew
        {
            //- Reference to the finite volume mesh this zone is part of
            const fvMesh& mesh_;

        public:

            iNew(const fvMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<porosityZone> operator()(Istream& is) const
            {
                word name(is);
                dictionary dict(is);

                return autoPtr<porosityZone>(new porosityZone(name, mesh_, dict));
            }
        };


    //- Destructor
    virtual ~porosityZone()
    {}


    // Member Functions

        // Access

            //- cellZone name
            const word& zoneName() const
            {
                return name_;
            }

            //- Return mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- cellZone number
            label zoneId() const
            {
                return cellZoneID_;
            }

            //- dictionary values used for the porosityZone
            const dictionary& dict() const
            {
                return dict_;
            }

            //- Return coordinate system
            const coordinateSystem& coordSys() const
            {
                return coordSys_;
            }

            //- Return origin
            const point& origin() const
            {
                return coordSys_.origin();
            }

            //- Return axis
#if EXTBRANCH==1
            vector axis() const
            {
                return coordSys_.axis();
            }
#elif OFPLUSBRANCH==1
            // NOTHING TO BE ADDED
#else
    #if OFVERSION<230
            vector axis() const
            {
                return coordSys_.axis();
            }
    #endif
#endif

            //- Return porosity
            scalar porosity() const
            {
                return porosity_;
            }

            //- Edit access to porosity
            scalar& porosity()
            {
                return porosity_;
            }

            //- Add the local porosities into the porosity field
            void porosity( volScalarField & ) const;

        //- Modify time derivative elements according to porosity
        template<class Type>
        void modifyDdt(fvMatrix<Type>&) const;

        //- Add the viscous and inertial resistance force contribution
        //  to the momentum equation
        void addResistance(fvVectorMatrix& UEqn) const;

        //- Add the viscous and inertial resistance force contribution
        //  to the tensorial diagonal.
        //  Optionally correct the processor BCs of AU.
        void addResistance
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU,
            bool correctAUprocBC = true
        ) const;

        //- Write the porosityZone dictionary
        virtual void writeDict(Ostream&, bool subDict = true) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const porosityZone&) = delete;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const porosityZone&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "porosityZoneTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
