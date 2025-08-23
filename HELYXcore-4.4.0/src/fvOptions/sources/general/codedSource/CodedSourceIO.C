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
    (c) 2012-2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "CodedSource.H"
#include "primitives/strings/stringOps/stringOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fv::CodedSource<Type>::sourceFields
(
    wordList& fieldNames
)
{
    fieldNames = coeffs_.lookup<wordList>("fields");
}


template<class Type>
bool Foam::fv::CodedSource<Type>::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        // Backward compatibility
        name_ = coeffs_.lookupOrDefault<word>("redirectType", "name");

        // Code snippets
        {
            const entry& e = coeffs_.lookupEntry
            (
                "codeCorrect",
                false,
                false
            );
            codeCorrect_ = stringOps::trim(e.stream());
            stringOps::inplaceExpand(codeCorrect_, coeffs_);
            dynamicCodeContext::addLineDirective
            (
                codeCorrect_,
                e.startLineNumber(),
                coeffs_.name()
            );
        }

        {
            const entry& e = coeffs_.lookupEntry
            (
                "codeAddSup",
                false,
                false
            );
            codeAddSup_ = stringOps::trim(e.stream());
            stringOps::inplaceExpand(codeAddSup_, coeffs_);
            dynamicCodeContext::addLineDirective
            (
                codeAddSup_,
                e.startLineNumber(),
                coeffs_.name()
            );
        }

        {
            const entry& e = coeffs_.lookupEntry
            (
                "codeSetValue",
                false,
                false
            );
            codeSetValue_ = stringOps::trim(e.stream());
            stringOps::inplaceExpand(codeSetValue_, coeffs_);
            dynamicCodeContext::addLineDirective
            (
                codeSetValue_,
                e.startLineNumber(),
                coeffs_.name()
            );
        }

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
