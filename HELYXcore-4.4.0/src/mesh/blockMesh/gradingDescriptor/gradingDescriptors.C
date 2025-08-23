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
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "gradingDescriptor/gradingDescriptors.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gradingDescriptors::gradingDescriptors()
:
    List<gradingDescriptor>(1, gradingDescriptor())
{}


Foam::gradingDescriptors::gradingDescriptors(const gradingDescriptor& gd)
:
    List<gradingDescriptor>(1, gd)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::gradingDescriptors Foam::gradingDescriptors::inv() const
{
    gradingDescriptors ret(*this);

    forAll(ret, i)
    {
        ret[i] = operator[](ret.size() - i - 1).inv();
    }

    return ret;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, gradingDescriptors& gds)
{
    // Examine next token
    token t(is);

    if (t.isNumber())
    {
        gds = gradingDescriptors(gradingDescriptor(t.number()));
    }
    else
    {
        is.putBack(t);

        // Read the list for gradingDescriptors
        is >> static_cast<List<gradingDescriptor>& >(gds);

        // Check state of Istream
        is.check(FUNCTION_NAME);

        // Normalize the blockFractions and nDivFractions
        // of the list of gradingDescriptors

        scalar sumBlockFraction = 0;
        scalar sumNDivFraction = 0;

        forAll(gds, i)
        {
            sumBlockFraction += gds[i].blockFraction_;
            sumNDivFraction += gds[i].nDivFraction_;
        }

        forAll(gds, i)
        {
            gds[i].blockFraction_ /= sumBlockFraction;
            gds[i].nDivFraction_ /= sumNDivFraction;
        }
    }

    return is;
}


// ************************************************************************* //
