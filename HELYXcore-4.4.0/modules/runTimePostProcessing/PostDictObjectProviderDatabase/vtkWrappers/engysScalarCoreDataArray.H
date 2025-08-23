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
    (c) 2025 Engys Ltd.

Class
    engysScalarCoreDataArray

Description
    Base wrapper for arrays belonging to the openFOAM core

SourceFiles
    engysScalarCoreDataArray.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "vtkWrappers/engysFoamCoreDataArray.H"
#include "vtkFloatArray.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
                   Class engysScalarCoreDataArray Declaration
\*---------------------------------------------------------------------------*/

class engysScalarCoreDataArray : public engysFoamCoreDataArray<Foam::scalar, Foam::scalar>
{
    typedef engysFoamCoreDataArray<Foam::scalar, Foam::scalar> engysFoamCoreDataArrayType;

public:
    vtkAbstractTypeMacroWithNewInstanceType(engysScalarCoreDataArray, engysFoamCoreDataArrayType, vtkDataArray, typeid(SelfType).name());

    static engysScalarCoreDataArray* New();

    engysScalarCoreDataArray() = default;

    engysScalarCoreDataArray(const engysScalarCoreDataArray&) = delete;
    void operator=(const engysScalarCoreDataArray&) = delete;

protected:
    VTK_NEWINSTANCE
    vtkObjectBase* NewInstanceInternal() const override { return vtkFloatArray::New(); }
};


// ************************************************************************* //
