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
    (c) 2019-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "SolverOption.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::SolverOption<Type>::SolverOption
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    solverOption(name, modelType, dict, obr),
    solverObject_(name, obr, dict)
{
    if (debug)
    {
        Info<< "call to SolverOption Constructor" << endl;
    }
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::SolverOption<Type>::~SolverOption()
{
    if (debug)
    {
        Info<< "call to SolverOption Destructor" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::fv::SolverOption<Type>::read(const dictionary& dict)
{
    if (dict.found("hookOp"))
    {
        word hookOp(dict.lookup("hookOp"));

        if
        (
            execHookTypeNames_[hookOp] == ehtCorrect
         || execHookTypeNames_[hookOp] == ehtConstrain
         || execHookTypeNames_[hookOp] == ehtOperator
        )
        {
            if (dict.found("hookField"))
            {
                word hookFieldName(dict.lookup("hookField"));
                bool found = false;
                forAll(fieldNames(), i)
                {
                    if (fieldNames()[i] == hookFieldName)
                    {
                        found = true;
                        break;
                    }
                }
                if (!found)
                {
                    fieldNames_.append(hookFieldName);
                    applied_.resize(fieldNames_.size(), false);
                }
            }
            else
            {
                FatalErrorInFunction
                    << "Mandatory option hookField not provided for hook "
                    << hookOp << " in fvOptions entry " << this->name() << nl
                    << exit(FatalError);
            }
        }
    }

    solverObject_.readBase(dict);
    return true;
}


template<class Type>
bool Foam::fv::SolverOption<Type>::initialise()
{
    if (debug)
    {
        Info<< "call to Foam::fv::solverOption::initialise()" << endl;
    }
    return solverObject_.initialiseBase();
}


template<class Type>
void Foam::fv::SolverOption<Type>::sourceFields
(
    wordList& fieldNames,
    wordList& regionNames
)
{
    solveList solveNames;
    SolveTable<solveList> sourceDependencies;
    solverObject_.getSourceGraphBase(solveNames, sourceDependencies);

    for (const solveID& solveName : solveNames)
    {
        fieldNames.append(solveName.first());
        regionNames.append(solveName.second());
    }
}


template<class Type>
void Foam::fv::SolverOption<Type>::boundarySourceFieldsAndPatches
(
    HashTable<labelList>& fieldPatchIDs
)
{
    HashTable<wordList> boundarySourceDependencies;
    solverObject_.getBoundarySourceGraph
    (
        fieldPatchIDs, boundarySourceDependencies
    );
}


template<class Type>
void Foam::fv::SolverOption<Type>::addSourceDependencies
(
    SolveTable<solveList>& dependencies
)
{
    solveList solveNames;
    SolveTable<solveList> sourceDependencies;
    solverObject_.getSourceGraphBase(solveNames, sourceDependencies);

    HashTable<labelList> boundarySourcePatchIDs;
    HashTable<wordList> boundarySourceDependencies;
    solverObject_.getBoundarySourceGraph
    (
        boundarySourcePatchIDs, boundarySourceDependencies
    );

    forAllConstIters(sourceDependencies, iter)
    {
        dependencies.insert(iter.key(), solveList());
        dependencies[iter.key()].append(iter());
    }

    word thisRegionName(obr_.name() == word::null ? mesh_.name() : obr_.name());
    forAllConstIters(boundarySourceDependencies, iter)
    {
        solveID key(iter.key(), thisRegionName);
        dependencies.insert(key, solveList());
        solveList& deps = dependencies[key];
        for (const word& dep : iter())
        {
            deps.append(solveID(dep, thisRegionName));
        }
    }
}

template<class Type>
void Foam::fv::SolverOption<Type>::correct()
{
    if (debug)
    {
        Info<< "call to Foam::fv::solverOption::correct()" << endl;
    }
    solverObject_.solve();
}


// ************************************************************************* //
