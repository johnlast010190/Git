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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "BlockKqRWallFunctionFvPatchField.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockKqRWallFunctionFvPatchField<Type>::BlockKqRWallFunctionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    kqRWallFunctionFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::BlockKqRWallFunctionFvPatchField<Type>::BlockKqRWallFunctionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    kqRWallFunctionFvPatchField<Type>(p, iF, dict)
{}


template<class Type>
Foam::BlockKqRWallFunctionFvPatchField<Type>::BlockKqRWallFunctionFvPatchField
(
    const BlockKqRWallFunctionFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    kqRWallFunctionFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::BlockKqRWallFunctionFvPatchField<Type>::BlockKqRWallFunctionFvPatchField
(
    const BlockKqRWallFunctionFvPatchField& tkqrwfpf
)
:
    kqRWallFunctionFvPatchField<Type>(tkqrwfpf)
{}


template<class Type>
Foam::BlockKqRWallFunctionFvPatchField<Type>::BlockKqRWallFunctionFvPatchField
(
    const BlockKqRWallFunctionFvPatchField& tkqrwfpf,
    const DimensionedField<Type, volMesh>& iF
)
:
    kqRWallFunctionFvPatchField<Type>(tkqrwfpf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
