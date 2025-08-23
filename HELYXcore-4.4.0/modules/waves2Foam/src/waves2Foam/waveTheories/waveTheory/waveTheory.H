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

Class
    Foam::waveTheories::waveTheory

Description
    An abstract base class for wave theories.

    A description of the general wave framework is given in

    @article { jacobsenFuhrmanFredsoe2011,
        Author = {Jacobsen, N G and Fuhrman, D R and Freds\o{}e, J},
        title = {{A Wave Generation Toolbox for the Open-Source CFD Library: OpenFoam\textregistered{}}},
        Journal = {{Int. J. for Numer. Meth. Fluids}},
        Year = {2012},
        Volume = {70},
        Number = {9},
        Pages = {1073-1088},
        DOI = {{10.1002/fld.2726}},
    }

SourceFiles
    waveTheory.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef waveTheory_H
#define waveTheory_H

#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "include/crossVersionCompatibility.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class waveTheory
:
    public IOdictionary
{
protected:

    // Protected data

        scalar seaLevel_;
        scalar rhoWater_;
        vector g_;
        vector direction_;

        const dictionary coeffDict_;

        scalar PI_;

        //- Wind vector
        vector wind_;


    // Protected member functions

        virtual scalar factor(const scalar&) const = 0;

        inline scalar returnZ(const point& x) const
        {
            return (-(direction_ & x) - seaLevel_);
        }

        void checkWaveDirection(const vector&) const;

        void checkWaveDirection(const vectorField&) const;


public:

    //- Runtime type information
    TypeName("waveTheory");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            waveTheory,
            dictionary,
            (
                const word& subDictName, const fvMesh& mesh_
            ),
            (subDictName, mesh_)
        );


    // Constructors

        //- Construct from components
        waveTheory
        (
            const word& type,
            const fvMesh& mesh_
        );

        //- Disallow default bitwise copy construction
        waveTheory(const waveTheory&) = delete;


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<waveTheory> New
        (
            const word& subDictName,
            const fvMesh& mesh_
        );

        virtual void printCoeffs() = 0;


    //- Destructor
    virtual ~waveTheory();


    // Member Functions

        //-
        virtual Switch conflictTSoftInitialise() const = 0;

        const vector windVelocity( const scalar& time ) const
        {
            return factor( time )*wind_;
        };

        // Access - single value
        virtual scalar eta
        (
            const point&,
            const scalar&
        ) const = 0;

        virtual scalar ddxPd
        (
            const point&,
            const scalar&,
            const vector&
        ) const = 0;

        virtual vector U
        (
            const point&,
            const scalar&
        ) const = 0;

        virtual scalar p
        (
            const point&,
            const scalar&
        ) const { return 0.0; };

        // Access - field
        virtual scalarField eta
        (
            const pointField&,
            const scalar&
        ) const;

        virtual scalarField ddxPd
        (
            const pointField&,
            const scalar&,
            const vectorField&
        ) const;

        virtual vectorField U
        (
            const pointField&,
            const scalar&
        ) const;

        //- Return the unit vertical direction in the direction of g
        const inline vector& returnDir() const
        {
            return direction_;
        }

        //- Return the sea level
        const inline scalar& returnSeaLevel() const
        {
            return seaLevel_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const waveTheory&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
