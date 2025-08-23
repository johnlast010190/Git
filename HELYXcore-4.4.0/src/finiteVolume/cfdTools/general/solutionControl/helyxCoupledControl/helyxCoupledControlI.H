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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2017 OpenCFD Ltd
    (c) 2020 Engys Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::helyxCoupledControl::nCorr() const
{
    return nCorr_;
}


inline Foam::label Foam::helyxCoupledControl::corr() const
{
    return corr_;
}


inline Foam::label& Foam::helyxCoupledControl::corr()
{
    return corr_;
}


inline Foam::label Foam::helyxCoupledControl::stabilityMode() const
{
    return stabilityMode_;
}


inline Foam::Switch Foam::helyxCoupledControl::choiCorrection() const
{
    return choiCorrection_;
}


inline bool Foam::helyxCoupledControl::correct()
{
    setFirstIterFlag();

    corr_++;

    if (debug)
    {
        Info<< algorithmName_ << " correct: corr = " << corr_ << endl;
    }

    if (corr_ <= nCorr_)
    {
        return true;
    }
    else
    {
        corr_ = 0;

        setFirstIterFlag();

        return false;
    }
}


inline bool Foam::helyxCoupledControl::storeInitialResiduals() const
{
    // Start from second PIMPLE iteration
    return corr_ == 2;
}


inline bool Foam::helyxCoupledControl::firstIter() const
{
    return corr_ == 1;
}


inline bool Foam::helyxCoupledControl::finalIter() const
{
    return converged_ || (corr_ == nCorr_);
}


inline bool Foam::helyxCoupledControl::solveFlow() const
{
    return solveFlow_;
}


inline bool Foam::helyxCoupledControl::turbCorr()
{
    if (turbOnFinalIterOnly_)
    {
        if (finalIter())
        {
            setFirstIterFlag(true, true);
        }
    }
    else
    {
        if (firstIter())
        {
            setFirstIterFlag(true, true);
        }
        else
        {
            setFirstIterFlag();
        }
    }

    return !turbOnFinalIterOnly_ || finalIter();
}


// ************************************************************************* //
