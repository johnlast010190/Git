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
    (c) 2025 Engys Ltd.

Application
    chemkinToFoam

Description
    Converts CHEMKINIII thermodynamics and reaction data files into
    OpenFOAM format.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "chemkinReader/chemkinReader.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "db/IOstreams/StringStreams/OStringStream.H"
#include "db/IOstreams/StringStreams/IStringStream.H"

#include "cfdTools/general/include/fvCFD.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Increase the precision of the output for JANAF coefficients
    Ostream::defaultPrecision(10);

    argList::addNote
    (
        "Converts CHEMKINIII thermodynamics and reaction data files into\n"
        "OpenFOAM format."
    );
    argList::noParallel();
    argList::noFunctionObjects();
    argList::validArgs.append("CHEMKINFile");
    argList::validArgs.append("CHEMKINThermodynamicsFile");
    argList::validArgs.append("CHEMKINTransport");
    argList::validArgs.append("FOAMChemistryFile");
    argList::validArgs.append("FOAMThermodynamicsFile");

    argList::addBoolOption
    (
        "newFormat",
        "read Chemkin thermo file in new format"
    );

    argList::addBoolOption
    (
        "legacyThermo",
        "write file in old thermodynamicProperties format"
    );

    argList args(argc, argv);

    #include "include/createTime.H"

    const word instance = runTime.findInstance
    (
        polyMesh::meshSubDir,
        "points",
        IOobject::READ_IF_PRESENT
    );

    const bool newFormat = args.optionFound("newFormat");
    const bool legacyThermo = args.optionFound("legacyThermo");

    chemkinReader cr(runTime, args[1], args[3], args[2], newFormat);

    const HashPtrTable<dictionary>& speciesThermo = cr.speciesThermo();

    dictionary thermoDict;
    thermoDict.add("species", cr.species());

    // Add the species thermo formatted entries
    {
        OStringStream os;
        forAllConstIter(HashPtrTable<dictionary>, speciesThermo, iter)
        {
            const dictionary& specieDict = iter();
            thermoDict.add(iter.key(), specieDict);
        }

        // Support of the new mateiral library introduced however it still
        // relies on the old thermo properties and just translates the strings
        // to new format. This will have to be properly converted before
        // we remove old thermo.
        if (legacyThermo)
        {
            const wordList toRemove
            ({
                "equationOfState",
                "thermodynamics",
                "kappaModel",
                "muModel",
                "molWeight"
            });
            forAllConstIter(HashPtrTable<dictionary>, speciesThermo, iter)
            {
                const word& specieName = iter()->name();
                dictionary& speciesDict = thermoDict.subDict(specieName);
                const scalar molWeight =
                    speciesDict.lookup<scalar>("molWeight");
                forAll(toRemove, i)
                {
                    speciesDict.remove(toRemove[i]);
                }
                dictionary specie("specie");
                specie.add("molWeight", molWeight);
                speciesDict.add("specie", specie);

                speciesDict.changeKeyword("muModelCoeffs", "transport");
                speciesDict.changeKeyword
                (
                    "thermodynamicsCoeffs",
                    "thermodynamics"
                );
            }
        }
        string speciesProperties = os.str();
        thermoDict.merge(dictionary(IStringStream(speciesProperties)()));
    }

    // Temporary hack to splice the specie composition data into the thermo file
    // pending complete integration into the thermodynamics structure

    // Add elements
    forAllConstIter(HashPtrTable<dictionary>, cr.speciesThermo(), iter)
    {
        const word specieName(iter.key());

        dictionary elementsDict("elements");
        forAll(cr.specieComposition()[specieName], ei)
        {
            elementsDict.add
            (
                cr.specieComposition()[specieName][ei].name(),
                cr.specieComposition()[specieName][ei].nAtoms()
            );
        }

        thermoDict.subDict(specieName).add("elements", elementsDict);
    }

    thermoDict.write(OFstream(args[5])(), false);

    cr.reactions().write(OFstream(args[4])());

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
