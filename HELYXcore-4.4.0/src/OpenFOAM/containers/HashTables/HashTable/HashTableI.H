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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "db/error/error.H"


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline Foam::label
Foam::HashTable<T, Key, Hash>::hashKeyIndex(const Key& key, label logicalSize) const
{
    // size is power of two - this is the modulus
    return Hash()(key) & (logicalSize - 1);
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline Foam::label Foam::HashTable<T, Key, Hash>::capacity() const
{
    return table_.size();
}


template<class T, class Key, class Hash>
inline Foam::label Foam::HashTable<T, Key, Hash>::size() const
{
    return nElmts_;
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::empty() const
{
    return !nElmts_;
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::insert
(
    const Key& key,
    const T& obj
)
{
    return this->set(key, obj, true);
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::set
(
    const Key& key,
    const T& obj
)
{
    return this->set(key, obj, false);
}


template<class T, class Key, class Hash>
inline const T& Foam::HashTable<T, Key, Hash>::lookup
(
    const Key& key,
    const T& deflt
) const
{
    const_iterator iter = this->find(key);
    return iter.found() ? iter.object() : deflt;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline T& Foam::HashTable<T, Key, Hash>::operator[](const Key& key)
{
    iterator iter = this->find(key);

    if (!iter.found())
    {
        FatalErrorInFunction
            << "Key \"" << key << "\" not found in table.  Valid entries: "
            << toc()
            << abort(FatalError);
    }

    return iter.object();
}


template<class T, class Key, class Hash>
inline const T& Foam::HashTable<T, Key, Hash>::operator[](const Key& key) const
{
    const_iterator iter = this->find(key);

    if (!iter.found())
    {
        FatalErrorInFunction
            << "Key \"" << key << "\" not found in table.  Valid entries: "
            << toc()
            << abort(FatalError);
    }

    return iter.object();
}


template<class T, class Key, class Hash>
inline T& Foam::HashTable<T, Key, Hash>::operator()(const Key& key)
{
    iterator iter = this->find(key);

    if (iter.found())
    {
        return iter.object();
    }

    this->insert(key, T());
    return find(key).object();
}


template<class T, class Key, class Hash>
inline T& Foam::HashTable<T, Key, Hash>::operator()
(
    const Key& key,
    const T& deflt
)
{
    iterator iter = this->find(key);

    if (iter.found())
    {
        return iter.object();
    }

    this->insert(key, deflt);
    return find(key).object();
}


template<class T, class Key, class Hash>
inline const T& Foam::HashTable<T, Key, Hash>::operator()
(
    const Key& key,
    const T& deflt
) const
{
    return this->lookup(key, deflt);
}


// * * * * * * * * * * * * * * * iterator base * * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline Foam::HashTable<T, Key, Hash>::iterator_base::iterator_base()
:
    entryPtr_(nullptr),
    hashTable_(nullptr),
    hashIndex_(0)
{}


template<class T, class Key, class Hash>
inline Foam::HashTable<T, Key, Hash>::iterator_base::iterator_base
(
    const table_type* hashTbl,
    const entry_type* elmt,
    const label hashIndex
)
:
    entryPtr_(const_cast<entry_type*>(elmt)),
    hashTable_(const_cast<table_type*>(hashTbl)),
    hashIndex_(hashIndex)
{}


template<class T, class Key, class Hash>
inline Foam::HashTable<T, Key, Hash>::iterator_base::iterator_base
(
    const table_type* hashTbl
)
:
    entryPtr_(nullptr),
    hashTable_(const_cast<table_type*>(hashTbl)),
    hashIndex_(0)
{
    if (hashTable_ && hashTable_->nElmts_)
    {
        // find first non-nullptr table entry
        while
        (
            !(entryPtr_ = hashTable_->table_[hashIndex_])
         && ++hashIndex_ < hashTable_->table_.size()
        )
        {}

        if (hashIndex_ >= hashTable_->table_.size())
        {
            // make into an end iterator
            entryPtr_ = nullptr;
            hashIndex_ = 0;
        }
    }
}


template<class T, class Key, class Hash>
inline void
Foam::HashTable<T, Key, Hash>::iterator_base::increment()
{
    // A negative index is a special value from erase
    if (hashIndex_ < 0)
    {
        // the markPos='-curPos-1', but we wish to continue at 'curPos-1'
        // thus use '-(markPos+1) -1'
        hashIndex_ = -(hashIndex_+1) - 1;
    }
    else if (entryPtr_)
    {
        if (entryPtr_->next_)
        {
            // Move to next element on the SLList
            entryPtr_ = entryPtr_->next_;
            return;
        }
    }
    // else
    // {
    //     // if we reach here (entryPtr_ is nullptr) it is already at the end()
    //     // we should probably stop
    // }


    // Step to the next table entry
    while
    (
        ++hashIndex_ < hashTable_->table_.size()
     && !(entryPtr_ = hashTable_->table_[hashIndex_])
    )
    {}

    if (hashIndex_ >= hashTable_->table_.size())
    {
        // make into an end iterator
        entryPtr_ = nullptr;
        hashIndex_ = 0;
    }
}


template<class T, class Key, class Hash>
inline bool
Foam::HashTable<T, Key, Hash>::iterator_base::found() const
{
    return entryPtr_;
}


template<class T, class Key, class Hash>
inline const Key& Foam::HashTable<T, Key, Hash>::iterator_base::key() const
{
    return entryPtr_->key_;
}


template<class T, class Key, class Hash>
inline T& Foam::HashTable<T, Key, Hash>::iterator_base::element() const
{
    if constexpr (HasData) {
        return entryPtr_->obj_;
    } else {
        return Foam::theNil;
    }
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::iterator_base::operator==
(
    const iterator_base& iter
) const
{
    return entryPtr_ == iter.entryPtr_;
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::iterator_base::operator!=
(
    const iterator_base& iter
) const
{
    return entryPtr_ != iter.entryPtr_;
}


// * * * * * * * * * * * * * * key iterator base * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
template<class WrappedIterator>
inline Foam::HashTable<T, Key, Hash>::key_iterator_base<WrappedIterator>
::key_iterator_base
(
    const WrappedIterator& iter
)
:
    WrappedIterator(iter)
{}


template<class T, class Key, class Hash>
template<class WrappedIterator>
inline const Key&
Foam::HashTable<T, Key, Hash>::key_iterator_base<WrappedIterator>
::operator*() const
{
    return this->key();
}


template<class T, class Key, class Hash>
template<class WrappedIterator>
inline const Key&
Foam::HashTable<T, Key, Hash>::key_iterator_base<WrappedIterator>
::operator()() const
{
    return this->key();
}


template<class T, class Key, class Hash>
template<class WrappedIterator>
inline typename Foam::HashTable<T, Key, Hash>::template key_iterator_base<WrappedIterator>&
Foam::HashTable<T, Key, Hash>::key_iterator_base<WrappedIterator>
::operator++()
{
    this->increment();
    return *this;
}


template<class T, class Key, class Hash>
template<class WrappedIterator>
inline typename Foam::HashTable<T, Key, Hash>::template key_iterator_base<WrappedIterator>
Foam::HashTable<T, Key, Hash>::key_iterator_base<WrappedIterator>
::operator++(int)
{
    key_iterator_base old = *this;
    this->increment();
    return old;
}

// ************************************************************************* //
