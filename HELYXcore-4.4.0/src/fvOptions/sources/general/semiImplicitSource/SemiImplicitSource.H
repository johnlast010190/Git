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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2017 Engys Ltd.

Class
    Foam::fv::SemiImplicitSource

Group
    grpFvOptionsSources

Description
    Time dependent semi-implicit source, described using an input dictionary.
    The injection rate coefficients are specified as separate files of Su-Sp
    coefficients, i.e.
        \f[
            S(x) = S_u + S_p x
        \f]
    where
    \vartable
        S(x)    | net source for field 'x'
        S_u     | explicit source contribution
        S_p     | linearised implicit contribution
    \endvartable

    Example of the source specification:

    \verbatim
    <Type>SemiImplicitSourceCoeffs
    {
        selectionMode   cellZone;
        cellZone        heater;

        volumeMode      absolute; // specific

        injectionRateSu //always constant <Type> or time vs <Type> table
        {
            k           table (( 0 0.005 ) (3.14159 0.01) (13.02 0.002));
            epsilon     constant 0.001;
            U           table (( 0 (1 1 1)) (1.22 (1.01 1.022 1.11)) );
            h           constant 101;
        }

        injectionRateSp //always constant <scalar> or time vs <scalars> table
        {
            k           constant 245;
            epsilon     table (( 0 0.005 ) (63.22 0.005 ) (65.11 0.000244 ));
            U           table (( 0 0.01 ) (6.12 0.0025 ) (11.1 0.244 ));
        }

        // alternative (old-style) interface
        injectionRateSuSp
        {
            k           (30.7 0);
            epsilon     (1.5  0);
        }
    }
    \endverbatim

    Notes
    - duplicate field names in the same list will give an error
    - if a field is only present in one list the other source
      for the field will be zero i.e. in above case h Sp source is zero

    Valid options for the \c volumeMode entry include:
    - absolute: values are given as \<quantity\>
    - specific: values are given as \<quantity\>/m3

SeeAlso
    Foam::fvOption

SourceFiles
    SemiImplicitSource.C

\*---------------------------------------------------------------------------*/

#ifndef SemiImplicitSource_H
#define SemiImplicitSource_H

#include "primitives/Tuple2/Tuple2.H"
#include "cellSetOption/cellSetOption.H"
#include "interpolations/interpolationTable/interpolationTable.H"
#include "containers/Lists/PtrList/PtrList.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "coordinate/systems/coordinateSystem.H"
#include "referenceFrames/coordinateFrame.H"
#include "coordinate/systems/cartesianCS.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

// Forward declaration of classes

template<class Type>
class SemiImplicitSource;


// Forward declaration of friend functions

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const SemiImplicitSource<Type>&
);


/*---------------------------------------------------------------------------*\
                     Class SemiImplicitSource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class SemiImplicitSource
:
    public cellSetOption
{
public:

    //- Private functions

        //- Source matrix
        template<class RhoFieldType>
        tmp<fvMatrix<Type>> sourceMatrix
        (
            const dimensionSet&,
            const VolField<Type>&,
            const RhoFieldType&,
            const label fieldi
        ) const;

        void createActivationField() const;

        void createSigmoidActivationField() const;

    // Public data

        //- Enumeration for volume types
        enum volumeModeType
        {
            vmAbsolute,
            vmSpecific
        };

        //- Word list of volume mode type names
        static const wordList volumeModeTypeNames_;


protected:

    // Protected data

        //- Volume mode
        volumeModeType volumeMode_;

        //- Volume normalisation
        mutable scalar VDash_;

        //- The time series being used for each variable
        PtrList<Function1<Type>> SuTimeSeriesList_;

        //- The time series being used for each variable
        PtrList<Function1<scalar>> SpTimeSeriesList_;

        bool useActivationField_;

        Switch useSigmoid_;

        scalar sigmoidFactor_;

        //- field generated by the Activation Function
        mutable autoPtr<volScalarField> activationField_;

        //- Activation Function
        autoPtr<Function1<vector>> activation_;

        FixedList<bool, 3> absoluteDistance_;

        //- Local coordinate System
        autoPtr<coordinateSystem> coorSysPtr_;

        //- Reference Frame Ptr
        coordinateFrame* coorFramePtr_;

        //- Apply density scaling y/n
        Switch densityScaling_;


    // Protected functions

        //- Helper function to convert from a word to a volumeModeType
        volumeModeType wordToVolumeModeType(const word& vtName) const;

        //- Helper function to convert from a volumeModeType to a word
        word volumeModeTypeToWord(const volumeModeType& vtType) const;

        //- find a word in a wordlist
        bool wordListFind(const wordList& wl, const word& w);

        //- Set the field time series
        void setTimeSeries(const dictionary& dict, wordList& fieldNames);

        scalar sigmoid(const scalar value, const scalar factor) const;

        //- Create coordinate system
        void createCS(const dictionary& dict);

        //- Local coordinate system
        inline const coordinateSystem& csys() const
        {
            if (coorFramePtr_)
            {
                return coorFramePtr_->coorSys();
            }
            return *coorSysPtr_;
        };

        //- Conversion from old interface (injectionRateSuSp)
        void convertOldInterface(const dictionary& dict, wordList& fieldNames);


public:

    //- Runtime type information
    TypeName("SemiImplicitSource");


    // Constructors

        //- Construct from components
        SemiImplicitSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    // Member Functions

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames);


        // Access

            //- Return const access to the volume mode
            inline const volumeModeType& volumeMode() const;

        // Edit

            //- Return access to the volume mode
            inline volumeModeType& volumeMode();

        // Evaluation

            //- Add explicit contribution to equation
            virtual void addSup
            (
                fvMatrix<Type>& eqn,
                const label fieldI
            );

            //- Add explicit contribution to compressible equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<Type>& eqn,
                const label fieldI
            );

            //- Add implicit contribution to block Momentum equation
            virtual void addSup
            (
                fvBlockMatrix<Type>& eqn,
                const label fieldi
            );

            //- Add implicit contribution to block Momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvBlockMatrix<Type>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SemiImplicitSource.C"
    #include "SemiImplicitSourceIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SemiImplicitSourceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
