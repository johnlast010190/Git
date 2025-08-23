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


SourceFiles
    thermalRunawaySource.C

\*---------------------------------------------------------------------------*/

#ifndef thermalRunawaySource_H
#define thermalRunawaySource_H

#include "cellSetOption/cellSetOption.H"
#include "db/functionObjects/stateFunctionObject/stateFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                         A templated CodedFvOption
\*---------------------------------------------------------------------------*/

class thermalRunawaySource
:
    public cellSetOption,
    public functionObjects::stateFunctionObject
{
public:

    //- Enumeration for volume types
    enum volumeModeType
    {
        vmAbsolute,
        vmSpecific
    };

    //- Word list of volume mode type names
    static const NamedEnum<volumeModeType, 2> volumeModeTypeNames_;

private:

    //- Volume mode
    volumeModeType volumeMode_;

    //- Volume normalisation
    scalar VDash_;

    //- Thermal runaway state
    bool state_;

    //- Time that thermal runaway starts
    scalar TR_time_;

    //- Critical temperature that initiates thermal runaway
    scalar Tcritical_;

    //- Heat Release Rate
    autoPtr<Function1<scalar>> HRR_;

public:

    //- Runtime type information
    TypeName("thermalRunawaySource");

    // Constructors

        //- Construct from patch and internal field
        thermalRunawaySource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        thermalRunawaySource(const thermalRunawaySource&) = delete;

        //- Destructor
        virtual ~thermalRunawaySource();


    // Member functions

        //- Return equations for which we provide a source
        void sourceFields(wordList& fieldNames) override;

        //- Explicit and implicit matrix contributions
        virtual void addSup
        (
            fvMatrix<scalar>& eqn,
            const label fieldi
        ) override;

        //- Explicit and implicit matrix contributions
        //- for compressible equations
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        ) override;

        //- Execute, does nothing
        bool execute() override;

        //- Write, does nothing
        bool write() override;

    // I-O

        //- Read dictionary
        virtual bool read(const dictionary& dict) override;


    // Member Operators

    //- Disallow default bitwise assignment
    void operator=(const thermalRunawaySource&) = delete;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


