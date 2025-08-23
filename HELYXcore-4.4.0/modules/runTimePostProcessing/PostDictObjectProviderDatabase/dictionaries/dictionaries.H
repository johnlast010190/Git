/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.0.1
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
    (c) 2015-2019 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::Dictionaries

Description
    Helper class to find the dictionaries used by the module

SourceFiles
    Dictionaries.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "db/Time/Time.H"
#include "db/dictionary/dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

#define GET_OPTIONAL_DICTIONARY(parent, key) (parent).isDict(key) ? (parent).subDict(key) : dictionary::null

/*---------------------------------------------------------------------------*\
                      Class Dictionaries Declaration
\*---------------------------------------------------------------------------*/

class Dictionaries
{
public:
    Dictionaries(
        const Time &runTime,
        const dictionary &postDict
    );

    [[nodiscard]] const dictionary& getRtppDict() const;
    [[nodiscard]] const dictionary& getHexMeshDict() const;
    [[nodiscard]] const dictionary& getBlockMeshDict() const;
    [[nodiscard]] const dictionary& getControlDict() const;
    [[nodiscard]] const dictionary& getFunctionObjectsDict() const;
    [[nodiscard]] const dictionary& getMeshObjectsDict() const;
    [[nodiscard]] const dictionary& getGeometryDict() const;

    ~Dictionaries();

private:
    dictionary hexMeshDict_;
    dictionary geometryDict_;
    dictionary blockMeshDict_;
    dictionary controlDict_;
    dictionary meshObjectsDict_;
    dictionary rtppDict_;

    static dictionary getHexMeshDict(const Time &runTime);
    static dictionary getGeometryDict(const Time &runTime);
    static dictionary getBlockMeshDict(const Time &runTime);
    static dictionary getControlDict(const Time &runTime);
    static dictionary getMeshObjectsDict(const Time &runTime);
    static dictionary getRtppDict(const Time &runTime, const dictionary& postDict);

    static dictionary getRequiredDictionary(const Time &runTime, const char *name);
    static dictionary getOptionalDictionary(const Time &runTime, const char *name);
    static dictionary getDictionary(const Time &runTime, const fileName &path);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
