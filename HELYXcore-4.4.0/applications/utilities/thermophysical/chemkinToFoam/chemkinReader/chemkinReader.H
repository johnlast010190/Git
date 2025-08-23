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
    (c) 2011-2023 OpenFOAM Foundation
    (c) 2024-2025 Engys Ltd.

Class
    Foam::chemkinReader

Group
    grpMulticomponentThermophysicalChemistryReaders

Description
    Foam::chemkinReader

SourceFiles
    chemkinReader.C
    chemkinLexer.C

\*---------------------------------------------------------------------------*/

#ifndef chemkinReader_H
#define chemkinReader_H

#include "primitives/strings/fileName/fileName.H"
#include "db/typeInfo/typeInfo.H"
#include "primitives/bools/Switch/Switch.H"
#include "containers/HashTables/HashPtrTable/HashPtrTable.H"
#include "reaction/Reactions/ReactionList/ReactionList.H"
#include "containers/Lists/DynamicList/DynamicList.H"
#include "primitives/ints/lists/labelList.H"
#include "speciesTable/speciesTable.H"
#include "specieElement/specieElement.H"
#include "atomicWeights/atomicWeights.H"
#include "include/dummyThermo.H"

#include <FlexLexer.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef HashTable<List<specieElement>> speciesCompositionTable;

/*---------------------------------------------------------------------------*\
                       Class chemkinReader Declaration
\*---------------------------------------------------------------------------*/

class chemkinReader
:
    public yyFlexLexer
{

public:

    // Public data types

        enum phase { solid, liquid, gas };

        //- Reference to object registry
        const objectRegistry& obr_;


private:

    // Private data


        static int yyBufSize;
        label lineNo_;

        //- Table of reaction type keywords
        HashTable<int> reactionKeywordTable_;

        //- Currently supported reaction types
        enum reactionKeyword
        {
            thirdBodyReactionType,
            unimolecularFallOffReactionType,
            chemicallyActivatedBimolecularReactionType,
            TroeReactionType,
            SRIReactionType,
            LandauTellerReactionType,
            reverseLandauTellerReactionType,
            JanevReactionType,
            powerSeriesReactionRateType,
            radiationActivatedReactionType,
            speciesTempReactionType,
            energyLossReactionType,
            plasmaMomentumTransfer,
            collisionCrossSection,
            nonEquilibriumReversibleReactionType,
            duplicateReactionType,
            speciesOrderForward,
            speciesOrderReverse,
            UnitsOfReaction,
            end
        };

        enum reactionType
        {
            irreversible,
            reversible,
            nonEquilibriumReversible,
            unknownReactionType
        };

        static const char* reactionTypeNames[4];

        enum reactionRateType
        {
            Arrhenius,
            thirdBodyArrhenius,
            unimolecularFallOff,
            chemicallyActivatedBimolecular,
            LandauTeller,
            Janev,
            powerSeries,
            unknownReactionRateType
        };

        static const char* reactionRateTypeNames[8];

        enum fallOffFunctionType
        {
            Lindemann,
            Troe,
            SRI,
            unknownFallOffFunctionType
        };

        static const char* fallOffFunctionNames[4];


        void initReactionKeywordTable();


        //- List of elements
        DynamicList<word> elementNames_;

        //- Element indices
        HashTable<label> elementIndices_;

        //- Isotope molecular weights
        HashTable<scalar> isotopeAtomicWts_;

        //- List of species
        DynamicList<word> specieNames_;

        //- Specie indices
        HashTable<label> specieIndices_;

        //- Table of species
        speciesTable speciesTable_;

        //- Specie phase
        HashTable<phase> speciePhase_;

        //- Table of the thermodynamic data given in the CHEMKIN file
        HashPtrTable<dictionary> speciesThermo_;

        //- List of the thermodynamic data given in the CHEMKIN file
        mutable PtrList<dummyThermo> speciesThermoList_;

        //- Table of species composition
        speciesCompositionTable speciesComposition_;

        //- List of the reactions
        ReactionList<dummyThermo> reactions_;

        //- Transport properties dictionary
        dictionary transportDict_;

        //- Flag to indicate that file is in new format
        Switch newFormat_;

        //- Tolerance for element imbalance in a reaction
        scalar imbalanceTol_;


    // Private Member Functions

        //- Flex lexer to read the CHEMKIN III file
        int lex();

        inline scalar stringToScalar(const string& s)
        {
            string& str = const_cast<string&>(s);
            str.replaceAll(" ", "");
            str.replaceAll("D", "e");
            str.replaceAll("d", "e");
            return atof(str.c_str());
        }

        inline scalar stringToScalar(const char* cstr)
        {
            return stringToScalar(string(cstr));
        }

        inline void correctElementName(word& elementName)
        {
            if (elementName.size() == 2)
            {
                elementName[1] = tolower(elementName[1]);
            }
            else if (elementName[0] == 'E')
            {
                elementName = "e";
            }
        }

        scalar molecularWeight
        (
            const List<specieElement>& specieComposition
        ) const;

        void finishElements(labelList& currentAtoms);

        void checkCoeffs
        (
            const scalarList& reactionCoeffs,
            const char* reationRateName,
            const label nCoeffs
        ) const;

        template<class ReactionRateType>
        void addReactionType
        (
            const reactionType rType,
            DynamicList<specieCoeffs>& lhs,
            DynamicList<specieCoeffs>& rhs,
            const ReactionRateType& rr
        );

        template<template<class, class> class PressureDependencyType>
        void addPressureDependentReaction
        (
            const reactionType rType,
            const fallOffFunctionType fofType,
            DynamicList<specieCoeffs>& lhs,
            DynamicList<specieCoeffs>& rhs,
            const scalarList& thirdBodyEfficiencies,
            const scalarList& k0Coeffs,
            const scalarList& kInfCoeffs,
            const HashTable<scalarList>& reactionCoeffsTable,
            const scalar Afactor0,
            const scalar AfactorInf,
            const scalar RR
        );

        void addReaction
        (
            DynamicList<specieCoeffs>& lhs,
            DynamicList<specieCoeffs>& rhs,
            const scalarList& thirdBodyEfficiencies,
            const reactionType rType,
            const reactionRateType rrType,
            const fallOffFunctionType fofType,
            const scalarList& ArrheniusReactionCoeffs,
            HashTable<scalarList>& reactionCoeffsTable,
            const scalar RR
        );

        // Read the CHEMKIN files
        void read
        (
            const fileName& CHEMKINFileName,
            const fileName& thermoFileName,
            const fileName& transportFileName
        );


public:

    // Constructors

        //- Construct from CHEMKIN III file name
        chemkinReader
        (
            const objectRegistry& obr,
            const fileName& CHEMKINFileName,
            const fileName& transportFileName,
            const fileName& thermoFileName = fileName::null,
            const bool newFormat = false
        );

        //- Disallow default bitwise copy construction
        chemkinReader(const chemkinReader&) = delete;


    //- Destructor
    virtual ~chemkinReader()
    {}


    // Member functions

        //- List of elements
        const wordList& elementNames() const
        {
            return elementNames_;
        }

        //- Element indices
        const HashTable<label>& elementIndices() const
        {
            return elementIndices_;
        }

        //- Isotope molecular weights
        const HashTable<scalar>& isotopeAtomicWts() const
        {
            return isotopeAtomicWts_;
        }

        //- Table of species
        const speciesTable& species() const
        {
            return speciesTable_;
        }

        //- Table of species composition
        const speciesCompositionTable& specieComposition() const
        {
            return speciesComposition_;
        }

        //- Specie phase
        const HashTable<phase>& speciePhase() const
        {
            return speciePhase_;
        }

        //- Table of the thermodynamic data given in the CHEMKIN file
        const HashPtrTable<dictionary>& speciesThermo() const
        {
            return speciesThermo_;
        }

        //- List of the thermodynamic data given in the CHEMKIN file
        const PtrList<dummyThermo>& speciesThermoList() const
        {
            speciesThermoList_.resize(speciesTable_.size());

            forAll(speciesTable_, speciei)
            {
                if (speciesThermoList_.set(speciei))
                {
                    continue;
                }

                HashPtrTable<dictionary>::const_iterator iter =
                    speciesThermo_.find(speciesTable_[speciei]);

                speciesThermoList_.set
                (
                    speciei,
                    new dummyThermo(obr_, iter())
                );
            }

            return speciesThermoList_;
        }

        //- List of the reactions
        const ReactionList<dummyThermo>& reactions() const
        {
            return reactions_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const chemkinReader&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
