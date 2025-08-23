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
    (c) 2016-2021 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.
    (c) 2019-2025 Engys Ltd.

Class
    Foam::functionObjects::specieReactionRates

Group
    grpFieldFunctionObjects

Description
    Writes the domain averaged, min, and max reaction rates for each species
    for each reaction into the file \<timeDir\>/specieReactionRates.dat,
    as well as the domain averaged, min and max heat release rates
    for each reaction

    If optional writeRRFields keyword is specified, write the reaction rate
    fields for the species and reactions specified as a list of tuples, i.e.
    writeRRFields ( (SpeciesName1 ReactionNo1) (SpeciesName2 ReactionNo2) ... )

    If optional writeQdotFields keyword is specified, write the heat release
    rate fields for the reactions specified in the list, i.e.
    writeQdotFields ( ReactionNo1 ReactionNo2 ... )

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::volRegion
    Foam::functionObjects::writeFile

SourceFiles
    specieReactionRates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_specieReactionRates_H
#define functionObjects_specieReactionRates_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "functionObjects/volRegion/volRegion.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "primitives/Tuple2/Tuple2.H"
#include "basicChemistryModel/basicChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class specieReactionRates Declaration
\*---------------------------------------------------------------------------*/

class specieReactionRates
:
    public fvMeshFunctionObject,
    public volRegion,
    public writeFile
{
    // Private Member Data

        const basicChemistryModel& chemistryModel_;

        //- Whether to save the Reaction Rate field per species / reaction
        List<Tuple2<word, label>> writeRRFields_;

        //- Whether to save the total heat release field per reaction
        labelList writeQdotFields_;


    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("specieReactionRates");


    // Constructors

        //- Construct from Time and dictionary
        specieReactionRates
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        specieReactionRates(const specieReactionRates&) = delete;


    //- Destructor
    virtual ~specieReactionRates();


    // Member Functions

        //- Read the specieReactionRates data
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Write the specie reaction rates
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const specieReactionRates&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
