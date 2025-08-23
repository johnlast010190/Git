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
    (c) 2021-2024 Engys Ltd.

Class
    Foam::gibFvPatchFieldMapper

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef gibFvPatchFieldMapper_H
#define gibFvPatchFieldMapper_H

#include "AMIInterpolation/AMIInterpolation/AMIInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class gibFvPatchFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class gibFvPatchFieldMapper
{
    private:

        const AMIInterpolation& amiInterp_;

    public:

        gibFvPatchFieldMapper(const AMIInterpolation& amiInterp)
        :
            amiInterp_(amiInterp)
        {}

        template<class Type>
        void
        map(Field<Type>& f, const Field<Type>& defaultValues) const
        {
            f.transfer
            (
                amiInterp_.interpolateToTarget
                (
                    f,
                    dynamic_cast<const UList<Type>&>(defaultValues)
                ).ref()
            );
        }

        template<class Type>
        void map(Field<Type>& f, const Type& defaultValue) const
        {
            Field<Type>
                defaultValues(amiInterp_.tgtMagSf().size(), defaultValue);
            map(f, defaultValues);
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
