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
    (c) 2021 Engys Ltd

Class
    Foam::HashTable

Description
    An STL-conforming hash table.

Note
    Hashing index collisions are handled via chaining using a singly-linked
    list with the colliding entry being added to the head of the linked
    list. Thus copying the hash table (or indeed even resizing it) will
    often result in a different hash order. Use a sorted table-of-contents
    when the hash order is important.

    The end iterator of all hash-tables has a nullptr to the hash entry.
    Thus avoid separate allocation for each table and use a single one with
    a nullptr. The hash-table iterators always have an entry-pointer as the
    first member data, which allows reinterpret_cast from anything else with
    a nullptr as its first data member.
    The nullObject is such an item (with a nullptr data member).

SourceFiles
    HashTableI.H
    HashTable.C
    HashTableCoreI.H
    HashTableCore.C
    HashTableIO.C

\*---------------------------------------------------------------------------*/

#ifndef HashTable_H
#define HashTable_H

#include "primitives/ints/label/label.H"
#include "primitives/ints/uLabel/uLabel.H"
#include "primitives/strings/word/word.H"
#include "primitives/hashes/Hash/Hash.H"
#include "db/typeInfo/className.H"
#include "primitives/nullObject/nullObject.H"
#include "containers/Lists/List/ListDecl.H"
#include "primitives/nil/nil.H"

#include <initializer_list>
#include <iterator>
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T, unsigned Size> class FixedList;
template<class T, class Key, class Hash> class HashTable;

template<class T, class Key, class Hash>
Istream& operator>>(Istream& is, HashTable<T, Key, Hash>& L);

template<class T, class Key, class Hash>
Ostream& operator<<(Ostream& os, const HashTable<T, Key, Hash>& tbl);


template<typename T, bool IsConst>
using maybe_const = std::conditional_t<IsConst, const T, T>;

/*---------------------------------------------------------------------------*\
                        Class HashTableCore Declaration
\*---------------------------------------------------------------------------*/

//- Bits that are independent of the HashTable template parameters.
struct HashTableCore
{
    //- Maximum allowable internal table size. Approximately labelMax/4
    static const label maxTableSize;

    //- Return a canonical (power-of-two) of the requested size.
    static label canonicalSize(const label requested_size);

    //- Define template name and debug
    ClassName("HashTable");

    //- Factory method to create a non-const iterator begin
    template<class IteratorType, class TableType>
    inline static IteratorType iterator_begin(TableType& table);

    //- Factory method to create a const iterator begin
    template<class IteratorType, class TableType>
    inline static IteratorType iterator_begin(const TableType& table);

    //- Factory method to create a const iterator begin
    template<class IteratorType, class TableType>
    inline static IteratorType iterator_cbegin(const TableType& table);

    //- Factory class for creating a begin/end pair for any const iterator.
    // TODO: KILL IT WITH FIRE.
    template<class IteratorType, class TableType>
    class const_iterator_pair
    {
        label size_;
        IteratorType iter_;
        IteratorType end_;

    public:

        inline const_iterator_pair(const TableType& tbl);

        inline label size() const;
        inline bool empty() const;

        inline IteratorType begin() const;
        inline IteratorType cbegin() const;

        inline const IteratorType end() const;
        inline const IteratorType cend() const;
    };

};


/*---------------------------------------------------------------------------*\
                           Class HashTable Declaration
\*---------------------------------------------------------------------------*/

template<class T, class Key=word, class Hash=string::hash>
class HashTable
:
    public HashTableCore
{
public:

    //- The template instance used for this HashTable
    typedef HashTable<T, Key, Hash> this_type;


    // STL type definitions

        //- Type of keys that the HashTable uses.
        typedef Key key_type;

        //- Type of values that the HashTable contains.
        typedef T value_type;

        //- The type used for storing into value_type objects.
        //  This type is usually value_type&.
        typedef T* pointer;

        //- The type used for storing into value_type objects.
        //  This type is usually value_type&.
        typedef T& reference;

        //- The type used for reading from constant value_type objects.
        typedef const T* const_pointer;

        //- The type used for reading from constant value_type objects.
        typedef const T& const_reference;

        //- The type to represent the difference between two iterators
        typedef label difference_type;

        //- The type that can represent the size of a HashTable.
        typedef label size_type;

        //- Forward-declare iterators
        template<bool IsConst>
        class iterator_impl;
        using iterator = iterator_impl<false>;
        using const_iterator = iterator_impl<true>;

private:

    // Private data type for table entries

        constexpr static bool HasData = !std::is_same_v<T, nil>;

        struct hashedEntryWithoutValues;
        struct hashedEntryWithValue;
        using hashedEntry = std::conditional_t<!HasData, hashedEntryWithoutValues, hashedEntryWithValue>;

        //- Structure to hold a hashed entry, with a linked-list for collisions
        struct hashedEntryWithoutValues
        {
            //- The lookup key
            Key key_;

            //- Pointer to next hashedEntry in sub-list
            hashedEntry* next_;

            //- Construct from key, object, next pointer
            hashedEntryWithoutValues(const Key& key, const T&, hashedEntry* next):
                key_(key), next_(next) {}

            hashedEntryWithoutValues(const Key& key, hashedEntry* next):
                key_(key), next_(next) {}
        };

        struct hashedEntryWithValue : public hashedEntryWithoutValues {
            //- The data object
            T obj_;

            //- Construct from key, object, next pointer
            hashedEntryWithValue(const Key& key, const T& obj, hashedEntry* next):
                hashedEntryWithoutValues(key, next),
                obj_(obj) {}
        };

    // Private data: size of table, the table and current number of elements

        //- The current number of elements in table
        label nElmts_;

        //- The table of primary entries
        List<hashedEntry*> table_;


    // Private Member Functions

        //- Return the hash index of the Key within the current table size.
        //  No checks for zero-sized tables.
        label hashKeyIndex(const Key& key, label logicalSize) const;
        label hashKeyIndex(const Key& key) const {
            return hashKeyIndex(key, table_.size());
        }

        //- Assign a new hash-entry to a possibly already existing key.
        //  \return True if the new entry was set.
        bool set(const Key& key, const T& obj, const bool protect);


protected:

        //- Internally used base for iterator and const_iterator
        class iterator_base;

        //- Friendship with the iterator_base is required.
        friend class iterator_base;


    // Protected Member Functions

        //- Remove using begin/end iterators of listed keys
        template<class InputIter>
        inline label eraseMultiple
        (
            const InputIter begIter,
            const InputIter endIter
        );


public:

    // Constructors

        //- Construct given initial table size
        HashTable(const label size = 128);

        //- Construct from Istream
        HashTable(Istream& is, const label size = 128);

        //- Construct as copy
        HashTable(const HashTable<T, Key, Hash>& ht);

        //- Move constructor by transferring the parameter contents
        HashTable(HashTable<T, Key, Hash>&& ht);

        //- Construct from an initializer list
        HashTable(std::initializer_list<std::pair<Key, T>> lst);


    //- Destructor
    ~HashTable();


    // Member Functions

      // Access

        //- The size of the underlying table
        inline label capacity() const;

        //- Return number of elements in table
        inline label size() const;

        //- Return true if the hash table is empty
        inline bool empty() const;

        //- Return true if hashed entry is found in table
        bool found(const Key& key) const;

        //- Find and return an iterator set at the hashed entry
        //  If not found iterator = end()
        iterator find(const Key& key);

        //- Find and return an const_iterator set at the hashed entry
        //  If not found iterator = end()
        const_iterator find(const Key& key) const;

        //- Find and return an const_iterator set at the hashed entry
        //  If not found iterator = end()
        const_iterator cfind(const Key& key) const;

        //- Return hashed entry if it exists, or return the given default
        inline const T& lookup(const Key& key, const T& deflt) const;


      // Table of contents

        //- Return the table of contents
        List<Key> toc() const;

        //- Return the table of contents as a sorted list
        List<Key> sortedToc() const;

        //- Return the sorted table of contents with keys that satisfy
        //  the unary predicate, optionally with inverted logic.
        template<class UnaryPredicate>
        List<Key> tocKeys
        (
            const UnaryPredicate& pred,
            const bool invert = false
        ) const;

        //- Return the sorted table of contents with values that satisfy
        // the unary predicate, optionally with inverted logic.
        template<class UnaryPredicate>
        List<Key> tocValues
        (
            const UnaryPredicate& pred,
            const bool invert = false
        ) const;

        //- Return the sorted table of contents with keys/values that satisfy
        //  the binary predicate, optionally with inverted logic.
        template<class BinaryPredicate>
        List<Key> tocEntries
        (
            const BinaryPredicate& pred,
            const bool invert = false
        ) const;


      // Counting

        //- Count the number of keys that satisfy the unary predicate,
        //  optionally with inverted logic.
        template<class UnaryPredicate>
        label countKeys
        (
            const UnaryPredicate& pred,
            const bool invert = false
        ) const;

        //- Count the number of values that satisfy the unary predicate,
        //  optionally with inverted logic.
        template<class UnaryPredicate>
        label countValues
        (
            const UnaryPredicate& pred,
            const bool invert = false
        ) const;

        //- Count the number of entries that satisfy the binary predicate,
        //  optionally with inverted logic.
        template<class BinaryPredicate>
        label countEntries
        (
            const BinaryPredicate& pred,
            const bool invert = false
        ) const;


      // Edit

        //- Insert a new entry
        //  \return True if the entry inserted, which means that it did
        //  not previously exist in the table.
        inline bool insert(const Key& key, const T& obj);

        //- Assign a new entry, overwriting existing entries.
        //
        //  \return True, since it always overwrites any entries.
        inline bool set(const Key& key, const T& obj);

        //- Erase an entry specified by given iterator
        //  This invalidates the iterator until the next ++ operation.
        //
        //  Includes a safeguard against the end-iterator such that the
        //  following is safe:
        //  \code
        //      auto iter = table.find(unknownKey);
        //      table.erase(iter);
        //  \endcode
        //  which is what \code table.erase(unknownKey) \endcode does anyhow.
        //
        //  \return True if the corresponding entry existed and was removed
        bool erase(const iterator& iter);

        //- Erase an entry specified by the given key
        //  \return True if the entry existed and was removed
        bool erase(const Key& key);

        //- Remove table entries given by the listed keys
        //  \return The number of items removed
        label erase(const UList<Key>& keys);

        //- Remove table entries given by the listed keys
        //  \return The number of items removed
        template<unsigned Size>
        label erase(const FixedList<Key, Size>& keys);

        //- Remove table entries given by the listed keys
        //  \return The number of items removed
        label erase(std::initializer_list<Key> keys);

        //- Remove table entries given by keys of the other hash-table.
        //
        //  The other hash-table must have the same type of key, but the
        //  type of values held and the hashing function are arbitrary.
        //
        //  \return The number of items removed
        template<class AnyType, class AnyHash>
        label erase(const HashTable<AnyType, Key, AnyHash>& other);

        //- Retain table entries given by keys of the other hash-table.
        //
        //  The other hash-table must have the same type of key, but the
        //  type of values held and the hashing function are arbitrary.
        //
        //  \return The number of items changed (removed)
        template<class AnyType, class AnyHash>
        label retain(const HashTable<AnyType, Key, AnyHash>& other);

        //- Generalized means to filter table entries based on their keys.
        //  Keep (or optionally prune) entries with keys that satisfy
        //  the unary predicate, which has the following signature:
        //  \code
        //  bool operator()(const Key& k);
        //  \endcode
        //
        //  For example,
        //  \code
        //  wordRes goodFields = ...;
        //  allFieldNames.filterKeys
        //  (
        //      [&goodFields](const word& k){ return goodFields.match(k); }
        //  );
        //  \endcode
        //
        //  \return The number of items changed (removed)
        template<class UnaryPredicate>
        label filterKeys
        (
            const UnaryPredicate& pred,
            const bool pruning = false
        );

        //- Generalized means to filter table entries based on their values.
        //  Keep (or optionally prune) entries with values that satisfy
        //  the unary predicate, which has the following signature:
        //  \code
        //  bool operator()(const T& v);
        //  \endcode
        //
        //  \return The number of items changed (removed)
        template<class UnaryPredicate>
        label filterValues
        (
            const UnaryPredicate& pred,
            const bool pruning = false
        );

        //- Generalized means to filter table entries based on their key/value.
        //  Keep (or optionally prune) entries with keys/values that satisfy
        //  the binary predicate, which has the following signature:
        //  \code
        //  bool operator()(const Key& k, const T& v);
        //  \endcode
        //
        //  \return The number of items changed (removed)
        template<class BinaryPredicate>
        label filterEntries
        (
            const BinaryPredicate& pred,
            const bool pruning = false
        );


        //- Resize the hash table for efficiency
        void resize(const label sz);

        //- Clear all entries from table
        void clear();

        //- Clear the table entries and the table itself.
        //  Equivalent to clear() followed by resize(0)
        void clearStorage();

        //- Transfer the contents of the argument table into this table
        //  and annul the argument table.
        void transfer(HashTable<T, Key, Hash>& ht);


    // Member Operators

        //- Find and return a hashed entry. FatalError if it does not exist.
        inline T& operator[](const Key& key);

        //- Find and return a hashed entry. FatalError if it does not exist.
        inline const T& operator[](const Key& key) const;

        //- Return existing entry or create a new entry.
        //  A newly created entry is created as a nameless T() and is thus
        //  value-initialized. For primitives, this will be zero.
        inline T& operator()(const Key& key);

        //- Return existing entry or insert a new entry.
        inline T& operator()(const Key& key, const T& deflt);

        //- Return hashed entry if it exists, or return the given default
        inline const T& operator()(const Key& key, const T& deflt) const;

        //- Assignment operator
        void operator=(const HashTable<T, Key, Hash>& rhs);

        //- Assignment operator from an initializer list
        void operator=(std::initializer_list<std::pair<Key, T>> lst);

        //- Move assignment operator
        void operator=(HashTable<T, Key, Hash>&&);

        //- Equality. Hash tables are equal if the keys and values are equal.
        //  Independent of table storage size and table order.
        bool operator==(const HashTable<T, Key, Hash>& rhs) const;

        //- The opposite of the equality operation. Takes linear time.
        bool operator!=(const HashTable<T, Key, Hash>& rhs) const;


protected:

    // Iterators and helpers

        class iterator_base
        {
        public:
            // Public typedefs
            using table_type = this_type;
            using key_type   = this_type::key_type;
            using iterator_category = std::forward_iterator_tag;
            using difference_type  = this_type::difference_type;

        private:
            using entry_type = hashedEntry;

          // Private Data

            //- Currently selected entry.
            //  MUST be the first member for easy comparison between iterators
            //  and for reinterpret_cast from nullObject
            entry_type* entryPtr_;

            //- Pointer to the hash-table for which this is an iterator
            //  This allows use of the default bitwise copy/assignment
            table_type* hashTable_;

            //- Current hash index within the hash-table data.
            //  A signed value, since erase() uses a negative value to signal
            //  the erasure state.
            label hashIndex_;

        protected:

          // Protected Member Functions

            //- Increment to the next position
            inline void increment();

            //- The referenced object/value element
            inline T& element() const;

            //- Erase the entry at the current position
            bool erase();


        public:

          // Constructors

            //- Construct null (end iterator)
            inline iterator_base();

            //- Construct from begin of hash-table
            inline explicit iterator_base(const table_type* hashTbl);

            //- Construct from hash table, element and hash index
            inline iterator_base
            (
                const table_type* hashTbl,
                const entry_type* elmt,
                const label hashIndex
            );

          // Member functions/operators

            //- True if iterator points to an entry
            //  This can be used directly instead of comparing to end()
            inline bool found() const;

            //- Return the Key corresponding to the iterator
            inline const Key& key() const;

            //- Compare hash-entry element pointers
            inline bool operator==(const iterator_base& iter) const;
            inline bool operator!=(const iterator_base& iter) const;
        };

public:

        //- An iterator wrapper for returning a reference to the key
        template<class WrappedIterator>
        class key_iterator_base
        :
            public WrappedIterator
        {
        public:
            using value_type = this_type::key_type;
            using pointer    = const Key*;
            using reference  = const Key&;

            //- Implicit conversion
            inline key_iterator_base(const WrappedIterator& iter);

            //- Return the key
            inline reference operator*() const;
            inline reference operator()() const;

            inline key_iterator_base& operator++();
            inline key_iterator_base operator++(int);
        };

    // STL iterator

        template<bool IsConst>
        class iterator_impl : public iterator_base {
            friend class HashTable;  // Uses iterator::erase() method
            using entry_type = maybe_const<hashedEntry, IsConst>;

            using iterator_t = iterator_impl<IsConst>;

        public:
            // Public typedefs
            using table_type = maybe_const<this_type, IsConst>;
            using value_type = maybe_const<this_type::value_type, IsConst>;
            using pointer    = maybe_const<this_type::pointer, IsConst>;
            using reference  = maybe_const<this_type::reference, IsConst>;

          // Constructors

            //- Construct null (end iterator)
            iterator_impl(): iterator_base() {}

            //- Construct from begin of hash-table
            explicit iterator_impl(table_type* hashTbl):
                iterator_base(hashTbl) {}

            //- Construct from hash table, element and hash index
            //  Used by the hash-table find() method.
            iterator_impl
            (
                table_type* hashTbl,
                entry_type* elmt,
                const label hashIndex
            ): iterator_base(hashTbl, elmt, hashIndex) {}

            // Non-const iterators can implicitly convert to const ones.
            template<bool OtherIsConst>
            iterator_impl(const iterator_impl<OtherIsConst>& other):
                iterator_base(other) {
                static_assert(OtherIsConst == IsConst || IsConst, "Cannot convert const iterator to non-const iterator");
            }

          // Member functions/operators

            //- Return non-const access to referenced object
            reference object() const {
                return this->element();
            }

            //- Return non-const access to referenced object
            reference operator*() const {
                return this->object();
            }
            reference operator()() const {
                return this->object();
            }

            iterator_t& operator++() {
                this->increment();
                return *this;
            }
            iterator_t operator++(int) {
                iterator_t old = *this;
                this->increment();
                return old;
            }
        };

    //- Iterating over keys only

        //- Forward iterator returning the key
        using key_iterator = key_iterator_base<iterator>;

        //- Forward const iterator returning the key
        using const_key_iterator = key_iterator_base<const_iterator>;

        //- A const iterator begin/end pair for iterating over keys
        const_iterator_pair<const_key_iterator, this_type> keys() const
        {
            return
                const_iterator_pair<const_key_iterator,this_type>(*this);
        }


    // Iterator access

        //- Iterator set to the beginning of the HashTable
        iterator begin() {
            return iterator(this);
        }

        //- const_iterator set to the beginning of the HashTable
        const_iterator begin() const {
            return const_iterator(this);
        }

        //- const_iterator set to the beginning of the HashTable
        const_iterator cbegin() const {
            return const_iterator(this);
        }

        //- iterator to signal the end for any HashTable
        iterator end() {
            return iterator{};
        }

        //- const_iterator to signal the end for any HashTable
        const_iterator end() const {
            return const_iterator{};
        }

        //- const_iterator to signal the end for any HashTable
        const_iterator cend() const {
            return const_iterator{};
        }


    // Writing

        //- Print information
        Ostream& printInfo(Ostream& os) const;

        //- Write the unordered keys as a list, with line-breaks if list length
        //  exceeds shortListLen. Using '0' suppresses line-breaks entirely.
        Ostream& writeKeys(Ostream& os, const label shortListLen=0) const;


    // IOstream Operator

        friend Istream& operator>> <T, Key, Hash>
        (
            Istream& is,
            HashTable<T, Key, Hash>& L
        );

        friend Ostream& operator<< <T, Key, Hash>
        (
            Ostream& os,
            const HashTable<T, Key, Hash>& tbl
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "containers/HashTables/HashTable/HashTableCoreI.H"
#include "containers/HashTables/HashTable/HashTableI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "containers/HashTables/HashTable/HashTable.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
