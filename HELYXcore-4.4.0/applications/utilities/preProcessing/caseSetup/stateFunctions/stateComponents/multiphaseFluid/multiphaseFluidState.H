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
    (c) 2016-2024 Engys Ltd.

Class
    Foam::multiphaseFluidState

Description
    Abstract, multiphase state common functions

SourceFiles
    multiphaseFluidState.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseFluidState_H
#define multiphaseFluidState_H

#include "stateComponents/turbulenceModel/turbulenceModelState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class multiphaseFluidState Declaration
\*---------------------------------------------------------------------------*/

class multiphaseFluidState
:
    public turbulenceModelState
{
protected:

    // Protected Data

        wordList phases_;


    // Protected member functions

        //- Translate generic phase field gradient entry to phase specific entry
        void addAlphaTypeGrad(const word& phaseName);

        //- Binary data assembly, normally scalar
        template<class Type>
        HashTable<Type, interfacePair, interfacePair::hash>
        assembleBinaryData
        (
            word property,
            const Type& defaultValue,
            const List<dictionary>&
        );

        //- Binary data assembly, dictionaries
        HashTable<dictionary, interfacePair, interfacePair::hash>
        assembleBinaryDictionaries
        (
            word property,
            const dictionary& defaultValue,
            const List<dictionary>&
        );

        //- Combine defaults with user input for material properties
        List<dictionary> assembleMaterialProperties();


public:

    //- Runtime type information
    TypeName("multiphaseFluid");

    // Constructors

        //- Construct from components
        multiphaseFluidState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );

        //- Disallow default bitwise copy construction
        multiphaseFluidState(const multiphaseFluidState&) = delete;


    //- Destructor
    virtual ~multiphaseFluidState()
    {}


    // Member Functions

        //- Remove default alpha entry from fieldMaps
        virtual void correct();

        //- Final manipulations
        virtual void finalise();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiphaseFluidState&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "multiphaseFluidStateTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
