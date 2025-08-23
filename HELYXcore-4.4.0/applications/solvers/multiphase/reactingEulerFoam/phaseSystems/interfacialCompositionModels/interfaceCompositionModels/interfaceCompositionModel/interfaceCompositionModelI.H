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
    (c) 2019-2020 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::phasePair& Foam::interfaceCompositionModel::pair() const
{
    return pair_;
}


const Foam::hashedWordList& Foam::interfaceCompositionModel::species() const
{
    return species_;
}


const Foam::rhoMulticomponentThermo&
Foam::interfaceCompositionModel::thermo() const
{
    return thermo_;
}


const Foam::basicSpecieMixture&
Foam::interfaceCompositionModel::composition() const
{
    return thermo_.composition();
}


const Foam::rhoThermo& Foam::interfaceCompositionModel::otherThermo() const
{
    return otherThermo_;
}


bool Foam::interfaceCompositionModel::otherHasComposition() const
{
    return isA<fluidMulticomponentThermo>(otherThermo_);
}


const Foam::basicSpecieMixture&
Foam::interfaceCompositionModel::otherComposition() const
{
    return refCast<const rhoMulticomponentThermo>(otherThermo_).composition();
}


// ************************************************************************* //