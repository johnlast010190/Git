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
    (c) 2010-2024 Engys Ltd.

Class
    Foam::fvSchemes

Description
    Selector class for finite volume differencing schemes.
    fvMesh is derived from fvShemes so that all fields have access to the
    fvSchemes from the mesh reference they hold.

SourceFiles
    fvSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef fvSchemes_H
#define fvSchemes_H

#include "db/IOobjects/IOdictionary/IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fvSchemes Declaration
\*---------------------------------------------------------------------------*/

class fvSchemes
:
    public IOdictionary
{
    // Private data

        dictionary ddtSchemes_;
        ITstream defaultDdtScheme_;

        dictionary d2dt2Schemes_;
        ITstream defaultD2dt2Scheme_;

        dictionary interpolationSchemes_;
        ITstream defaultInterpolationScheme_;

        dictionary divSchemes_;
        ITstream defaultDivScheme_;

        dictionary gradSchemes_;
        ITstream defaultGradScheme_;

        dictionary snGradSchemes_;
        ITstream defaultSnGradScheme_;

        dictionary laplacianSchemes_;
        ITstream defaultLaplacianScheme_;

        mutable dictionary fluxRequired_;
        bool defaultFluxRequired_;

    //- Dictionary of direct set fvScheme parameters
    dictionary localSchemes_;

        //- Steady-state run indicator
        //  Set true if the default ddtScheme is steadyState
        bool steady_;


    // Private Member Functions

        //- Clear the dictionaries and streams before reading
        void clear();


public:
        //- Read settings from the dictionary
        void read(const dictionary&);


private:

        //- Refresh fvSchemes dependent member data
        //- and apply local modifications
        void update();


public:

    //- Debug switch
    static int debug;


    // Constructors

        //- Construct for objectRegistry
        fvSchemes(const objectRegistry& obr);

        //- Disallow default bitwise copy construction
        fvSchemes(const fvSchemes&) = delete;


    // Member Functions

        // Access

            const dictionary& dict() const;

            ITstream& ddtScheme(const word& name) const;

            ITstream& d2dt2Scheme(const word& name) const;

            ITstream& interpolationScheme(const word& name) const;

            ITstream& divScheme(const word& name) const;

            ITstream& gradScheme(const word& name) const;

            ITstream& snGradScheme(const word& name) const;

            ITstream& laplacianScheme(const word& name) const;

            const dictionary& ddtSchemes() const
            {
                return ddtSchemes_;
            }

            const dictionary& d2dt2Schemes() const
            {
                return d2dt2Schemes_;
            }

            const dictionary& interpolationSchemes() const
            {
                return interpolationSchemes_;
            }

            const dictionary& divSchemes() const
            {
                return divSchemes_;
            }

            const dictionary& gradSchemes() const
            {
                return gradSchemes_;
            }

            const dictionary& snGradSchemes() const
            {
                return snGradSchemes_;
            }

            const dictionary& laplacianSchemes() const
            {
                return laplacianSchemes_;
            }

            void setFluxRequired(const word& name) const;

            bool fluxRequired(const word& name) const;

            //- Copy input to localSchemes_
            // Call will force a reread of this and apply
            // localSchemes via a reverse merge operation to this
            // after every read
            void setLocalSchemeDict(const dictionary&);

            //- const acces to local dictionary
            const dictionary& localSchemeDict() const
            {
                return localSchemes_;
            }

            void resetReadOpt(const readOption& option);

            //- Return true if the default ddtScheme is steadyState
            bool steady() const
            {
                return steady_;
            }

            //- Return true if the default ddtScheme is not steadyState
            bool transient() const
            {
                return !steady_;
            }


        // Read

            //- Read the fvSchemes
            bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvSchemes&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
