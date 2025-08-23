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
    (c) 2013-2020 OpenFOAM Foundation

Class
    Foam::radiationModels::sootModels::mixtureFraction

Group
    grpRadiationSootSubModels

Description
    This soot model is purely an state model. The amount of soot produced is
    determined by a single step chemistry as :

        nuf Fuel + nuOx Ox = nuP P + nuSoot soot

    nuSoot is prescribed by the user.

    The single step chemistry used is read from the combustion.
    The soot is not considered into the thermodynamics of the system and it
    is not considered as an extra specie in the solver.

    The spacial distribution is given by the normalization of the first product
    on the rhs of the reaction by default or it can be added as input.

    For example in the radiationProperties dictionary set:

    sootModel mixtureFraction;

    mixtureFractionCoeffs
    {
        nuSoot              0.015;
        Wsoot               12;
        mappingField        P;
    }

SourceFiles
    mixtureFraction.C

\*---------------------------------------------------------------------------*/

#ifndef mixtureFraction_H
#define mixtureFraction_H

#include "interpolations/interpolationLookUpTable/interpolationLookUpTable.H"
#include "sootModels/sootModel/sootModel.H"
#include "containers/HashTables/HashTable/HashTable.H"
#include "fluidThermo/fluidThermo.H"
#include "mixtures/multicomponentMixture/multicomponentMixture.H"
#include "reaction/Reactions/Reaction/Reaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace sootModels
{

/*---------------------------------------------------------------------------*\
                    Class mixtureFraction Declaration
\*---------------------------------------------------------------------------*/
class mixtureFraction
:
    public sootModel
{
    // Private Data

         //- Soot mass fraction
        volScalarField soot_;

        //- Soot model dictionary
        dictionary coeffsDict_;

        //- Soot yield
        scalar nuSoot_;

        //- Soot molecular weight
        scalar Wsoot_;

        //- Maximum soot mass concentration at stoichiometric
        scalar sootMax_;

        //- Name of the field mapping the soot
        word mappingFieldName_;

        //- Maximum value of the map field
        scalar mapFieldMax_;


public:

    //- Runtime type information
    TypeName("mixtureFraction");


    // Constructors

        //- Construct from components
        mixtureFraction
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& modelType
        );


    //- Destructor
    virtual ~mixtureFraction();


    // Member Functions

        //- Return Ysoot
        const volScalarField& soot() const
        {
            return soot_;
        }

        //- Main update/correction routine
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sootModels
} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
