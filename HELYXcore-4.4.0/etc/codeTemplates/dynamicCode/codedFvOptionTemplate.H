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
    (c) 2022 Engys Ltd.

Description
    Template for use with dynamic code generation of a source.
    The hook functions take the following arguments:

    codeCorrect
    (
        VolField<Type>& fld
    )

    codeAddSup
    (
        fvMatrix<Type}>& eqn,
        const label fieldi
    )

    setValue
    (
        fvMatrix<Type}>& eqn,
        const label fieldi
    )

    where :
        fieldi is the index in the fields entry
        eqn is the fvMatrix

    energySource
    {
        type            scalarCodedSource;

        active          true;
        selectionMode   all;

        scalarCodedSourceCoeffs
        {
            fields  (h);
            name    sourceTime;

            codeInclude
            #{

            #};

            codeCorrect
            #{
                Pout<< "**codeCorrect**" << endl;
            #};

            codeAddSup
            #{
                const Time& time = mesh().time();
                const scalarField& V = mesh_.V();
                scalarField& heSource = eqn.source();
                heSource -= 0.1*sqr(time.value())*V;
            #};

            codeSetValue
            #{
                Pout<< "**codeSetValue**" << endl;
            #};

            // Dummy entry. Make dependent on above to trigger recompilation
            code
            #{
                $codeInclude
                $codeCorrect
                $codeAddSup
                $codeSetValue
            #};
        }

        sourceTimeCoeffs
        {
            // Dummy entry
        }
    }

SourceFiles
    codedFvOptionTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef codedFvOptionTemplate_H
#define codedFvOptionTemplate_H

#include "cellSetOption/cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
                         A templated CodedFvOption
\*---------------------------------------------------------------------------*/

class ${typeName}FvOption${SourceType}
:
    public cellSetOption
{
public:

    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;

    //- Runtime type information
    TypeName("${typeName}");


    // Constructors

        //- Construct from patch and internal field
        ${typeName}FvOption${SourceType}
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

    //- Destructor
    virtual ~${typeName}FvOption${SourceType}();


    // Member functions

            //- Correct field
            virtual void correct
            (
                VolField<${TemplateType}>&
            );

            //- Explicit and implicit matrix contributions
            virtual void addSup
            (
                fvMatrix<${TemplateType}>& eqn,
                const label fieldi
            );

            //- Explicit and implicit matrix contributions for compressible
            //  equations
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<${TemplateType}>& eqn,
                const label fieldi
            );

            //- Set value
            virtual void setValue
            (
                fvMatrix<${TemplateType}>& eqn,
                const label fieldi
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
