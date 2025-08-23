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
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::ColourMap

Description
    Data storage for colour map data

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "db/dictionary/dictionary.H"
#include "primitives/strings/word/word.H"
#include "primitives/Scalar/scalar/scalar.H"
#include "primitives/Vector/vector/vector.H"
#include "postDict/postDictKeys.H"
#include "types/colourSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{
// Structs to hold data
class ColourMap
{
// This is basically a data structure - no need to be fussy about access
public:

    word name;
    ColourSpace colourSpace;
    bool hsvWrap = false;
    List<scalar> points;
    List<Vector<scalar>> colours;

    // Construct from dictionary
    ColourMap() = default;
    explicit ColourMap(const dictionary& dict)
    {
        name = dict.lookup<word>(idKeys::NAME_KEY);
        colourSpace = dict.lookup<ColourSpace>(colourMapKeys::COLOR_SPACE_KEY);
        hsvWrap = dict.lookup<bool>(colourMapKeys::HSV_WRAP_KEY);
        points = dict.lookup<List<scalar>>(colourMapKeys::POINTS_KEY);
        colours = dict.lookup<List<Vector<scalar>>>(colourMapKeys::COLORS_KEY);
    }

    [[nodiscard]] scalar validPointForIndex(label index) const
    {
        if (hasValidPoints())
        {
            return points[index];
        }
        else
        {
            return static_cast<scalar>(index) / static_cast<scalar>(colours.size() - 1);
        }
    }

    [[nodiscard]] bool hasValidPoints() const
    {
        return points.size() == colours.size();
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
