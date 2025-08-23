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
    (c) 2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::linearEqn::linearEqn()
{}


inline Foam::linearEqn::linearEqn(const Foam::zero)
:
    VectorSpace<linearEqn, scalar, 2>(Foam::zero())
{}


inline Foam::linearEqn::linearEqn(const scalar a, const scalar b)
{
    this->v_[A] = a;
    this->v_[B] = b;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::linearEqn::a() const
{
    return this->v_[A];
}


inline Foam::scalar Foam::linearEqn::b() const
{
    return this->v_[B];
}


inline Foam::scalar& Foam::linearEqn::a()
{
    return this->v_[A];
}


inline Foam::scalar& Foam::linearEqn::b()
{
    return this->v_[B];
}


inline Foam::scalar Foam::linearEqn::value(const scalar x) const
{
    return x*a() + b();
}


inline Foam::scalar Foam::linearEqn::error(const scalar x) const
{
    return mag(SMALL*x*a());
}


inline Foam::Roots<1> Foam::linearEqn::roots() const
{
    /*

    This function solves a linear equation of the following form:

        a*x + b = 0
          x + B = 0

    */

    const scalar a = this->a();
    const scalar b = this->b();

    if (a == 0)
    {
        return Roots<1>(roots::nan, 0);
    }

    if (mag(b/VGREAT) >= mag(a))
    {
        return Roots<1>(sign(a) == sign(b) ? roots::negInf : roots::posInf, 0);
    }

    return Roots<1>(roots::real, - b/a);
}


// ************************************************************************* //
