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
    Foam::relaxationWeights::relaxationWeight

Description
    An abstract base class for generating different relaxation Weights through
    a dictionary lookup.

    See e.g.

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
    relaxationWeight.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef relaxationWeight_H
#define relaxationWeight_H

#include "fields/volFields/volFields.H"

#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "relaxationWeightCourantCorrection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationWeights
{

/*---------------------------------------------------------------------------*\
                      Class relaxationShape Declaration
\*---------------------------------------------------------------------------*/

class relaxationWeight
:
    public IOdictionary
{
protected:

    // Protected data

        const fvMesh& mesh_;

        const dictionary coeffDict_;

        scalar PI_;

        relaxationWeightCourantCorrection rwcc_;


    // Protected member functions

        virtual void computeWeights
        (
            const labelList&,
            const scalarField&,
            scalarField&
        ) = 0;


public:

    //- Runtime type information
    TypeName("relaxationWeight");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            relaxationWeight,
            dictionary,
            (
                const word& subDictName, const fvMesh& mesh_
            ),
            (subDictName, mesh_)
        );


    // Constructors

        //- Construct from components
        relaxationWeight
        (
            const word& type,
            const fvMesh& mesh_
        );

        //- Disallow default bitwise copy construction
        relaxationWeight(const relaxationWeight&) = delete;


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<relaxationWeight> New
        (
            const word& subDictName,
            const fvMesh& mesh_
        );


    //- Destructor
    virtual ~relaxationWeight();


    // Member Functions

        // Access
        virtual void weights
        (
            const labelList&,
            const scalarField&,
            scalarField&
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const relaxationWeight&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationWeights
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
