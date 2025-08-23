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
    (c) 2007 Takuya OSHIMA <oshima@eng.niigata-u.ac.jp>
    (c) 2024 Engys Ltd.

Description
    Convective outlet boundary condition.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/convectiveOutlet/convectiveOutletFvPatchField.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
convectiveOutletFvPatchField<Type>::convectiveOutletFvPatchField
(
    const fvPatch& p,
#if FIELD_IS_DIMENSIONED
    const DimensionedField<Type, volMesh>& iF
#else
    const Field<Type>& iF
#endif
)
:
    fixedGradientFvPatchField<Type>(p, iF),
    convectiveVelocity_(p.size(), 0.0),
    curTimeIndex_(0),
    snGradScheme_("normal"),
    ddtScheme_("CrankNicholson"),
    updateValue_(false),
    writeValue_(false),
    fieldPatchName_("::" + p.name()),
    gi0_(p.size(), pTraits<Type>::zero),
    gb0_(p.size(), pTraits<Type>::zero),
    pi0_(p.size(), pTraits<Type>::zero),
    pi00_(p.size(), pTraits<Type>::zero),
    pb0_(p.size(), pTraits<Type>::zero),
    pb00_(p.size(), pTraits<Type>::zero)
{}


template<class Type>
convectiveOutletFvPatchField<Type>::convectiveOutletFvPatchField
(
    const convectiveOutletFvPatchField<Type>& ptf,
    const fvPatch& p,
#if FIELD_IS_DIMENSIONED
    const DimensionedField<Type, volMesh>& iF,
#else
    const Field<Type>& iF,
#endif
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchField<Type>(ptf, p, iF, mapper),
    convectiveVelocity_(mapper(ptf.convectiveVelocity_)),
    curTimeIndex_(ptf.curTimeIndex_),
    snGradScheme_(ptf.snGradScheme_),
    ddtScheme_(ptf.ddtScheme_),
    updateValue_(ptf.updateValue_),
    writeValue_(ptf.writeValue_),
    fieldPatchName_(ptf.fieldPatchName_),
    gi0_(mapper(ptf.gi0_)),
    gb0_(mapper(ptf.gb0_)),
    pi0_(mapper(ptf.pi0_)),
    pi00_(mapper(ptf.pi00_)),
    pb0_(mapper(ptf.pb0_)),
    pb00_(mapper(ptf.pb00_))
{}


template<class Type>
convectiveOutletFvPatchField<Type>::convectiveOutletFvPatchField
(
    const fvPatch& p,
#if FIELD_IS_DIMENSIONED
    const DimensionedField<Type, volMesh>& iF,
#else
    const Field<Type>& iF,
#endif
    const dictionary& dict
)
:
    fixedGradientFvPatchField<Type>(p, iF),
    convectiveVelocity_("convectiveVelocity", dict, p.size()),
    curTimeIndex_(0),
    snGradScheme_("normal"),
    ddtScheme_("CrankNicholson"),
    updateValue_(false),
    writeValue_(false),
    fieldPatchName_(dict.name().name()),
    gi0_(p.size(), pTraits<Type>::zero),
    gb0_(p.size(), pTraits<Type>::zero),
    pi0_(p.size(), pTraits<Type>::zero),
    pi00_(p.size(), pTraits<Type>::zero),
    pb0_(p.size(), pTraits<Type>::zero),
    pb00_(p.size(), pTraits<Type>::zero)
{
    if (dict.found("snGradScheme"))
    {
        snGradScheme_ = dict.lookup<word>("snGradScheme");

        if
        (
            snGradScheme_ != "upwind"
         && snGradScheme_ != "predictorCorrector"
         && snGradScheme_ != "normal"
        )
        {
            FatalErrorIn
            (
                "convectiveOutletFvPatchField::convectiveOutletFvPatchField()"
            )   << "    Unsupported surface-normal differencing scheme : "
                << snGradScheme_ << exit(FatalError);
        }
    }

    if (dict.found("ddtScheme"))
    {
        ddtScheme_ = dict.lookup<word>("ddtScheme");

        if
        (
            ddtScheme_ != "Euler"
         && ddtScheme_ != "backward"
         && ddtScheme_ != "CrankNicholson"
        )
        {
            FatalErrorIn
            (
                "convectiveOutletFvPatchField::convectiveOutletFvPatchField()"
            )   << "    Unsupported temporal differencing scheme : "
                << ddtScheme_ << exit(FatalError);
        }
    }

    if (dict.found("updateValue"))
    {
        updateValue_ = dict.lookup<bool>("updateValue");
    }

    if (dict.found("writeValue"))
    {
        writeValue_ = dict.lookup<bool>("writeValue");
    }

    pi00_ = this->patchInternalField();
    pi0_ = this->patchInternalField();
    if (dict.found("gradient"))
    {
        this->gradient() = Field<Type>("gradient", dict, p.size());
        pb00_ = pi00_ + this->gradient()/this->patch().deltaCoeffs();
    }
    else
    {
        this->gradient() = Field<Type>(this->size(), pTraits<Type>::zero);
        pb00_ = pi00_;
    }
    pb0_ = pb00_;

    this->forceAssign(this->patchInternalField()
  + this->gradient()/this->patch().deltaCoeffs());
}


template<class Type>
convectiveOutletFvPatchField<Type>::convectiveOutletFvPatchField
(
    const convectiveOutletFvPatchField<Type>& coptf,
#if FIELD_IS_DIMENSIONED
    const DimensionedField<Type, volMesh>& iF
#else
    const Field<Type>& iF
#endif
)
:
    fixedGradientFvPatchField<Type>(coptf, iF),
    convectiveVelocity_(coptf.convectiveVelocity_),
    curTimeIndex_(coptf.curTimeIndex_),
    snGradScheme_(coptf.snGradScheme_),
    ddtScheme_(coptf.ddtScheme_),
    updateValue_(coptf.updateValue_),
    writeValue_(coptf.writeValue_),
    fieldPatchName_(coptf.fieldPatchName_),
    gi0_(coptf.gi0_),
    gb0_(coptf.gb0_),
    pi0_(coptf.pi0_),
    pi00_(coptf.pi00_),
    pb0_(coptf.pb0_),
    pb00_(coptf.pb00_)
{}

#if NEEDS_COPY_CTOR
template<class Type>
convectiveOutletFvPatchField<Type>::convectiveOutletFvPatchField
(
    const convectiveOutletFvPatchField<Type>& coptf
)
    :
    fixedGradientFvPatchField<Type>(coptf),
    convectiveVelocity_(coptf.convectiveVelocity_),
    curTimeIndex_(coptf.curTimeIndex_),
    snGradScheme_(coptf.snGradScheme_),
    ddtScheme_(coptf.ddtScheme_),
    updateValue_(coptf.updateValue_),
    writeValue_(coptf.writeValue_),
    fieldPatchName_(coptf.fieldPatchName_),
    gi0_(coptf.gi0_),
    gb0_(coptf.gb0_),
    pi0_(coptf.pi0_),
    pi00_(coptf.pi00_),
    pb0_(coptf.pb0_),
    pb00_(coptf.pb00_)
{}
#endif

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void convectiveOutletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const Time& ldb = this->db().time();

    if (ddtScheme_ == "Euler")
    {
        if (snGradScheme_ == "normal")
        {
            this->gradient() =
                (pb0_ - this->patchInternalField())
               /(convectiveVelocity_*ldb.deltaT().value()
              + 1.0/this->patch().deltaCoeffs());
        }
        else if (snGradScheme_ == "upwind")
        {
            this->gradient() =
                -(this->patchInternalField() - pi0_)
                /(convectiveVelocity_*ldb.deltaT().value());
        }
        else if (snGradScheme_ == "predictorCorrector")
        {
            // Predictor
            this->gradient() =
              - (this->patchInternalField() - pi0_)
               /(convectiveVelocity_*ldb.deltaT().value());

            // Corrector
            this->gradient() =
                0.5*
                (
                    this->gradient()
                  - (
                        (
                            this->patchInternalField()
                          + this->gradient()/this->patch().deltaCoeffs()
                        )
                      - pb0_
                    )
                   /(convectiveVelocity_*ldb.deltaT().value())
                );
        }
    }
    else if (ddtScheme_ == "backward")
    {
        if (snGradScheme_ == "normal")
        {
            this->gradient() =
                (4.0*pb0_ - pb00_ - 3.0*this->patchInternalField())
               /(
                    2.0*convectiveVelocity_*ldb.deltaT().value()
                  + 3.0/this->patch().deltaCoeffs()
                );
        }
        else if (snGradScheme_ == "upwind")
        {
            this->gradient() =
               -(3.0*this->patchInternalField() - 4.0*pi0_ + pi00_)
               /(2.0*convectiveVelocity_*ldb.deltaT().value());
        }
        else if (snGradScheme_ == "predictorCorrector")
        {
            // Predictor
            this->gradient() =
              - (3.0*this->patchInternalField() - 4.0*pi0_ + pi00_)
               /(2.0*convectiveVelocity_*ldb.deltaT().value());

            // Corrector
            this->gradient() =
                0.5*
                (
                    this->gradient()
                  - (
                        3.0
                       *(
                            this->patchInternalField()
                          + this->gradient()/this->patch().deltaCoeffs()
                        )
                        - 4.0*pb0_ + pb00_
                    )/(2.0*convectiveVelocity_*ldb.deltaT().value())
                );
        }
    }
    else if (ddtScheme_ == "CrankNicholson")
    {
        if (snGradScheme_ == "normal")
        {
            this->gradient() =
                (
                    2.0*(pb0_ - this->patchInternalField())
                  - convectiveVelocity_*ldb.deltaT().value()*gb0_
                )
               /(
                    convectiveVelocity_*ldb.deltaT().value()
                  + 2.0/this->patch().deltaCoeffs()
                );
        }
        else if (snGradScheme_ == "upwind")
        {
            this->gradient() =
                -2.0
               *(this->patchInternalField() - pi0_)
               /(convectiveVelocity_*ldb.deltaT().value())
              - gb0_;
        }
        else if (snGradScheme_ == "predictorCorrector")
        {
            // Predictor
            this->gradient() =
                -2.0
                *(this->patchInternalField() - pi0_)
                /(convectiveVelocity_*ldb.deltaT().value()) - gi0_;

            // Use time index to save oldTime gradient values
            if (curTimeIndex_ != ldb.timeIndex())
            {
                gi0_ = this->gradient();
            }

            // Corrector
            this->gradient() =
                0.5*
                (
                    this->gradient()
                  - 2.0*
                    (
                        (
                            this->patchInternalField()
                          + this->gradient()/this->patch().deltaCoeffs()
                        )
                      - pb0_
                    )
                   /(convectiveVelocity_*ldb.deltaT().value()) - gb0_
                );
        }
    }

    if (updateValue_)
    {
        this->forceAssign
        (
            this->patchInternalField()
          + this->gradient()/this->patch().deltaCoeffs()
        );
    }

    // Use time index to save oldTime patchField values
    if (curTimeIndex_ != ldb.timeIndex())
    {
        gb0_ = this->gradient();

        if (ddtScheme_ == "backward")
        {
            pi00_ = pi0_;
            pb00_ = pb0_;
        }
        pi0_ = this->patchInternalField();

        pb0_ =
            this->patchInternalField()
          + this->gradient()/this->patch().deltaCoeffs();

        curTimeIndex_ = ldb.timeIndex();
    }

    // Sets Updated to true
    fixedGradientFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void convectiveOutletFvPatchField<Type>::write(Ostream& os) const
{
    fixedGradientFvPatchField<Type>::write(os);
    os.writeEntry("snGradScheme", snGradScheme_);
    os.writeEntry("ddtScheme", ddtScheme_);
    os.writeEntry("updateValue", updateValue_);
    os.writeEntry("writeValue", writeValue_);
    convectiveVelocity_.writeEntry("convectiveVelocity", os);

    if (writeValue_)
    {
        this->writeEntry("value", os);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
