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
    (c) 2011-20121 OpenFOAM Foundation

Class
    Foam::polynomialTransport

Group
    grpSpecieTransport

Description
    Transport package using polynomial functions for \c mu and \c kappa.

Usage

    \table
        Property        | Description
        muCoeffs<8>     | Dynamic viscosity polynomial coefficients
        kappaCoeffs<8>  | Thermal conductivity polynomial coefficients
    \endtable

    Example of the specification of the transport properties:
    \verbatim
    transport
    {
        muCoeffs<8>     ( 1000 -0.05 0.003 0 0 0 0 0 );
        kappaCoeffs<8>  ( 2000 -0.15 0.023 0 0 0 0 0 );
    }
    \endverbatim

    The polynomial expressions are evaluated as so:

        \f[
            \mu    = 1000 - 0.05 T + 0.003 T^2
        \f]

        \f[
            \kappa = 2000 - 0.15 T + 0.023 T^2
        \f]

SourceFiles
    polynomialTransportI.H
    polynomialTransport.C

See also
    Foam::Polynomial

\*---------------------------------------------------------------------------*/

#ifndef polynomialTransport_H
#define polynomialTransport_H

#include "primitives/functions/Polynomial/Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo, int PolySize> class polynomialTransport;

template<class Thermo, int PolySize>
inline polynomialTransport<Thermo, PolySize> operator+
(
    const polynomialTransport<Thermo, PolySize>&,
    const polynomialTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
inline polynomialTransport<Thermo, PolySize> operator*
(
    const scalar,
    const polynomialTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
Ostream& operator<<
(
    Ostream&,
    const polynomialTransport<Thermo, PolySize>&
);

class objectRegistry;

/*---------------------------------------------------------------------------*\
                     Class polynomialTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, int PolySize=8>
class polynomialTransport
:
    public Thermo
{
    // Private data

        //- Dynamic viscosity polynomial coefficients [Pa.s/K^i]
        Polynomial<PolySize> muCoeffs_;

        //- Thermal conductivity polynomial coefficients [W/m/K/K^i]
        Polynomial<PolySize> kappaCoeffs_;


    // Private Member Functions

        //- Construct from components
        inline polynomialTransport
        (
            const Thermo& t,
            const Polynomial<PolySize>& muPoly,
            const Polynomial<PolySize>& kappaPoly
        );


public:

    // Constructors

        //- Copy constructor
        polynomialTransport(const polynomialTransport&) = default;

        //- Construct as named copy
        inline polynomialTransport(const word&, const polynomialTransport&);

        //- Construct from dictionary
        polynomialTransport(const objectRegistry& obr, const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<polynomialTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<polynomialTransport> New
        (
            const objectRegistry& obr,
            const dictionary& dict
        );

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

            autoPtr<polynomialTransport> operator()(const dictionary& dict) const
            {
                return New(mesh_, dict);
            }
        };

    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "polynomial<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline void operator=(const polynomialTransport&);

        inline void operator+=(const polynomialTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend polynomialTransport operator+ <Thermo, PolySize>
        (
            const polynomialTransport&,
            const polynomialTransport&
        );

        friend polynomialTransport operator* <Thermo, PolySize>
        (
            const scalar,
            const polynomialTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo, PolySize>
        (
            Ostream&,
            const polynomialTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "transport/polynomial/polynomialTransportI.H"

#ifdef NoRepository
    #include "transport/polynomial/polynomialTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
