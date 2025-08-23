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
    (c) 2019 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include <type_traits>
#include "vtk/output/foamVtkOutput.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::vtk::fileWriter::writeUniform
(
    const word& fieldName,
    const Type& val,
    const label nValues
)
{
    static_assert
    (
        (
            std::is_same<label, typename pTraits<Type>::cmptType>::value
         || std::is_floating_point<typename pTraits<Type>::cmptType>::value
        ),
        "Label and Floating-point vector space only"
    );

    const direction nCmpt(pTraits<Type>::nComponents);

    if (format_)
    {
        if (std::is_same<label, typename pTraits<Type>::cmptType>::value)
        {
            if (legacy())
            {
                legacy::intField<nCmpt>(format(), fieldName, nValues);
            }
            else
            {
                const uint64_t payLoad = vtk::sizeofData<label, nCmpt>(nValues);

                format().beginDataArray<label, nCmpt>(fieldName);
                format().writeSize(payLoad);
            }
        }
        else
        {
            if (legacy())
            {
                legacy::floatField<nCmpt>(format(), fieldName, nValues);
            }
            else
            {
                const uint64_t payLoad = vtk::sizeofData<float, nCmpt>(nValues);

                format().beginDataArray<float, nCmpt>(fieldName);
                format().writeSize(payLoad);
            }
        }
    }

    if (format_)
    {
        for (label i=0; i < nValues; ++i)
        {
            vtk::write(format(), val);
        }
    }

    if (format_)
    {
        format().flush();
        format().endDataArray();
    }
}


// ************************************************************************* //
