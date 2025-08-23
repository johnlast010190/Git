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
    (c) 2022 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "codedFvOptionTemplate.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "global/unitConversion/unitConversion.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"

//{{{ begin codeInclude
${codeInclude}
//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode
${localCode}
//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = ${SHA1sum}
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void ${typeName}_${SHA1sum}(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//makeRemovablePatchTypeField
//(
//    fvPatch${FieldType},
//    ${typeName}FvOption${SourceType}
//);
defineTypeNameAndDebug(${typeName}FvOption${SourceType}, 0);
addRemovableToRunTimeSelectionTable
(
    option,
    ${typeName}FvOption${SourceType},
    dictionary
);


const char* const ${typeName}FvOption${SourceType}::SHA1sum =
    "${SHA1sum}";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

${typeName}FvOption${SourceType}::
${typeName}FvOption${SourceType}
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr)
{
    if (${verbose:-false})
    {
        Info<<"construct ${typeName} sha1: ${SHA1sum}"
            " from components\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

${typeName}FvOption${SourceType}::
~${typeName}FvOption${SourceType}()
{
    if (${verbose:-false})
    {
        Info<<"destroy ${typeName} sha1: ${SHA1sum}\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ${typeName}FvOption${SourceType}::correct
(
    VolField<${TemplateType}>& fld
)
{
    if (${verbose:-false})
    {
        Info<<"${typeName}FvOption${SourceType}::correct()\n";
    }

//{{{ begin code
    ${codeCorrect}
//}}} end code
}


void ${typeName}FvOption${SourceType}::addSup
(
    fvMatrix<${TemplateType}>& eqn,
    const label fieldi
)
{
    if (${verbose:-false})
    {
        Info<<"${typeName}FvOption${SourceType}::addSup()\n";
    }

//{{{ begin code
    ${codeAddSup}
//}}} end code
}


void ${typeName}FvOption${SourceType}::addSup
(
    const volScalarField& rho,
    fvMatrix<${TemplateType}>& eqn,
    const label fieldi
)
{
    if (${verbose:-false})
    {
        Info<<"${typeName}FvOption${SourceType}::addSup()\n";
    }

//{{{ begin code
    ${codeAddSup}
//}}} end code
}


void ${typeName}FvOption${SourceType}::setValue
(
    fvMatrix<${TemplateType}>& eqn,
    const label fieldi
)
{
    if (${verbose:-false})
    {
        Info<<"${typeName}FvOption${SourceType}::setValue()\n";
    }

//{{{ begin code
    ${codeSetValue}
//}}} end code
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

} // End namespace fv
// ************************************************************************* //
