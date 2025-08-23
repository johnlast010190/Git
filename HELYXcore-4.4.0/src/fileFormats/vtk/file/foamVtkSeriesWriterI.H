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
    (c) 2018-2019 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::vtk::seriesWriter::empty() const
{
    return entries_.empty();
}


inline Foam::label Foam::vtk::seriesWriter::size() const
{
    return entries_.size();
}


inline void Foam::vtk::seriesWriter::clear()
{
    entries_.clear();
    existing_.clear();
}


inline bool Foam::vtk::seriesWriter::append(const fileNameInstant& inst)
{
    // Strip out path before saving
    return appendCheck(fileNameInstant(inst.value(), inst.name().name()));
}


inline bool Foam::vtk::seriesWriter::append(fileNameInstant&& inst)
{
    // Strip out path before saving
    inst.name().removePath();

    return appendCheck(inst);
}


inline bool Foam::vtk::seriesWriter::append
(
    scalar timeValue,
    const fileName& file
)
{
    // Strip out path before saving
    return appendCheck(fileNameInstant(timeValue, file.name()));
}


inline bool Foam::vtk::seriesWriter::append
(
    scalar timeValue,
    fileName&& file
)
{
    // Strip out path before saving
    file.removePath();

    return appendCheck(fileNameInstant(timeValue, std::move(file)));
}


inline void Foam::vtk::seriesWriter::print(Ostream& os) const
{
    seriesWriter::print(os, entries_);
}


inline void Foam::vtk::seriesWriter::write(const fileName& seriesName) const
{
    seriesWriter::write(seriesName, entries_);
}


// ************************************************************************* //
