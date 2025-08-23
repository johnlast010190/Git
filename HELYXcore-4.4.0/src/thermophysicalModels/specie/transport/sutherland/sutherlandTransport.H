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
    (c) 2011-2021 OpenFOAM Foundation

Class
    Foam::sutherlandTransport

Group
    grpSpecieTransport

Description
    Transport package using Sutherland's formula for viscosity:

    \verbatim
        mu = As*sqrt(T)/(1 + Ts/T)
    \endverbatim

    and the thermal conductivity (for gases) is obtained from:

    \verbatim
        kappa = mu*Cv*(1.32 + 1.77*R/Cv)
    \endverbatim

Usage
    \table
        Property     | Description
        As           | Sutherland constant kg/(ms K^1/2)
        Ts           | Sutherland temperature [K]
    \endtable

    Example specification of sutherlandTransport for air:

    \verbatim
    transport
    {
        As              1.458e-06;
        Ts              110.4;
    }
    \endverbatim

SourceFiles
    sutherlandTransportI.H
    sutherlandTransport.C

\*---------------------------------------------------------------------------*/

#ifndef sutherlandTransport_H
#define sutherlandTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class sutherlandTransport;

template<class Thermo>
inline sutherlandTransport<Thermo> operator+
(
    const sutherlandTransport<Thermo>&,
    const sutherlandTransport<Thermo>&
);

template<class Thermo>
inline sutherlandTransport<Thermo> operator*
(
    const scalar,
    const sutherlandTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const sutherlandTransport<Thermo>&
);

class fvMesh;

/*---------------------------------------------------------------------------*\
                    Class sutherlandTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class sutherlandTransport
:
    public Thermo
{
    // Private data

        // Sutherland's coefficients
        scalar As_, Ts_;


    // Private Member Functions

        //- Calculate the Sutherland coefficients
        //  given two viscosities and temperatures
        inline void calcCoeffs
        (
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
        );

        //- Read coefficient from dictionary
        scalar readCoeff(const word& coeffName, const dictionary& dict);


public:

    // Constructors

        //- Copy constructor
        sutherlandTransport(const sutherlandTransport&) = default;

        //- Construct from components
        inline sutherlandTransport
        (
            const Thermo& t,
            const scalar As,
            const scalar Ts
        );

        //- Construct from two viscosities
        inline sutherlandTransport
        (
            const Thermo& t,
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
        );

        //- Construct as named copy
        inline sutherlandTransport(const word&, const sutherlandTransport&);

        //- Construct from dictionary
        sutherlandTransport(const objectRegistry& obr, const dictionary& dict);

        //- Construct from base thermo and dictionary
        sutherlandTransport(const Thermo& t, const objectRegistry& mesh, const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<sutherlandTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<sutherlandTransport> New(const objectRegistry& mesh, const dictionary& dict);

        //- Return a pointer to a new phaseModel created on freestore
        //  from Istream
        class iNew
        {
            const objectRegistry& mesh_;

        public:

            iNew
            (
                const objectRegistry& mesh
            )
            :
                mesh_(mesh)
            {}

            autoPtr<sutherlandTransport> operator()(const dictionary& dict) const
            {
                return New(mesh_, dict);
            }
        };

    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "sutherland<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline void operator=(const sutherlandTransport&);

        inline void operator+=(const sutherlandTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend sutherlandTransport operator+ <Thermo>
        (
            const sutherlandTransport&,
            const sutherlandTransport&
        );

        friend sutherlandTransport operator* <Thermo>
        (
            const scalar,
            const sutherlandTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const sutherlandTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "transport/sutherland/sutherlandTransportI.H"

#ifdef NoRepository
    #include "transport/sutherland/sutherlandTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
