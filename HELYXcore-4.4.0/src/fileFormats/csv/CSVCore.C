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
    (c) 2011 OpenFOAM Foundation
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "csv/CSVCore.H"
#include "db/IOstreams/StringStreams/IStringStream.H"
#include "containers/Lists/DynamicList/DynamicList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::CSVCore::CSVCore()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::List<Foam::string> Foam::fileFormats::CSVCore::readLine
(
    const string& line,
    const char& separator,
    const Switch& mergeSeparators,
    const label& nEntries
)
{
    //- code copied from CSV.C and modified
    label n = 0;
    std::size_t pos = 0;

    DynamicList<string> splitted;

    if (mergeSeparators)
    {
        std::size_t nPos = 0;

        while ((pos != std::string::npos) && (n <= nEntries))
        {
            bool found = false;
            while (!found)
            {
                nPos = line.find(separator, pos);

                if ((nPos != std::string::npos) && (nPos - pos == 0))
                {
                    pos = nPos + 1;
                }
                else
                {
                    found = true;
                }
            }

            nPos = line.find(separator, pos);

            if (nPos == std::string::npos)
            {
                splitted.append(line.substr(pos));
                pos = nPos;
                n++;
            }
            else
            {
                splitted.append(line.substr(pos, nPos - pos));
                pos = nPos + 1;
                n++;
            }
        }
    }
    else
    {
        while ((pos != std::string::npos) && (n <= nEntries))
        {
            std::size_t nPos = line.find(separator, pos);

            if (nPos == std::string::npos)
            {
                splitted.append(line.substr(pos));
                pos = nPos;
                n++;
            }
            else
            {
                splitted.append(line.substr(pos, nPos - pos));
                pos = nPos + 1;
                n++;
            }
        }
    }
    return List<string>(splitted, true);
}


// ************************************************************************* //
