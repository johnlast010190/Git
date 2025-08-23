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

\*---------------------------------------------------------------------------*/

inline Foam::scalar Foam::polynomialFunction::value
(
    const scalar x
) const
{
    const scalarList& coeffs = *this;
    scalar val = coeffs[0];

    // avoid costly pow() in calculation
    scalar powX = x;
    for (label i=1; i<coeffs.size(); ++i)
    {
        val += coeffs[i]*powX;
        powX *= x;
    }

    if (logActive_)
    {
        val += this->logCoeff_*log(x);
    }

    return val;
}

inline Foam::scalar Foam::polynomialFunction::derivative
(
    const scalar x
) const
{
    const scalarList& coeffs = *this;
    scalar deriv = 0;
    label polySize=coeffs.size();
    if (polySize > 1)
    {
        // avoid costly pow() in calculation
        deriv += coeffs[1];

        scalar powX = 1;
        for (label i=2; i<polySize; ++i)
        {
            powX *= x;
            deriv += i*coeffs[i]*powX;
        }
    }

    if (logActive_)
    {
        deriv += logCoeff_/x;
    }

    return deriv;

}

