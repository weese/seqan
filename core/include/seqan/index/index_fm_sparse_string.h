// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#ifndef INDEX_FM_SPARSE_STRING_H_
#define INDEX_FM_SPARSE_STRING_H_

namespace seqan {

// ==========================================================================
// Forwards
// ==========================================================================

template <typename TString, typename TSpec>
struct SparseString;

// ==========================================================================
// Tags
// ==========================================================================

// ----------------------------------------------------------------------------
// Tag FibreValues
// ----------------------------------------------------------------------------

struct FibreValues_;
typedef Tag<FibreValues_>       const FibreValues;

// ----------------------------------------------------------------------------
// Tag FibreIndicators
// ----------------------------------------------------------------------------

struct FibreIndicators_;
typedef Tag<FibreIndicators_>   const FibreIndicators;

// ==========================================================================
// Metafunctions
// ==========================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TFibreValues, typename TSpec>
struct Value<SparseString<TFibreValues, TSpec> >
{
    typedef typename Value<TFibreValues>::Type Type;
};

template <typename TFibreValues, typename TSpec>
struct Value<SparseString<TFibreValues, TSpec> const> :
    Value<SparseString<TFibreValues, TSpec> > {};

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

template <typename TFibreValues, typename TSpec>
struct GetValue<SparseString<TFibreValues, TSpec> > :
    Value<SparseString<TFibreValues, TSpec> > {};

template <typename TFibreValues, typename TSpec>
struct GetValue<SparseString<TFibreValues, TSpec> const> :
    Value<SparseString<TFibreValues, TSpec> const> {};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TFibreValues, typename TSpec>
struct Reference<SparseString<TFibreValues, TSpec> >
{
    typedef typename Value<SparseString<TFibreValues, TSpec> >::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction DefaultValue
// ----------------------------------------------------------------------------

template <typename TSpec>
struct DefaultValue;

template <typename TSpec>
struct DefaultValue<TSpec const> :
    DefaultValue<TSpec> {};

template <typename TFibreValues, typename TSpec>
struct DefaultValue<SparseString<TFibreValues, TSpec> >
{
    typedef typename GetValue<SparseString<TFibreValues, TSpec> const>::Type Type;
    static const Type VALUE = -1;
};

// ----------------------------------------------------------------------------
// Metafunction Fibre
// ----------------------------------------------------------------------------

template <typename TFibreValues, typename TSpec>
struct Fibre<SparseString<TFibreValues, TSpec>, FibreValues>
{
    typedef TFibreValues Type;
};

template <typename TFibreValues, typename TSpec>
struct Fibre<SparseString<TFibreValues, TSpec>, FibreIndicators>
{
    typedef RankDictionary<TwoLevels<bool, TSpec> > Type;
};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

template <typename TFibreValues, typename TSpec>
struct Iterator<SparseString<TFibreValues, TSpec>, Standard>
{
    typedef Iter<SparseString<TFibreValues, TSpec>, PositionIterator> Type;
};

template <typename TFibreValues, typename TSpec>
struct Iterator<SparseString<TFibreValues, TSpec> const, Standard>
{
    typedef Iter<SparseString<TFibreValues, TSpec> const, PositionIterator> Type;
};

template <typename TFibreValues, typename TSpec>
struct Iterator<SparseString<TFibreValues, TSpec>, Rooted>:
    Iterator<SparseString<TFibreValues, TSpec>, Standard>{};

template <typename TFibreValues, typename TSpec>
struct Iterator<SparseString<TFibreValues, TSpec> const, Rooted>:
    Iterator<SparseString<TFibreValues, TSpec> const, Standard>{};

// ==========================================================================
// Classes
// ==========================================================================

// ----------------------------------------------------------------------------
// Class SparseString
// ----------------------------------------------------------------------------
// NOTE(esiragusa): Why is SparseString not a specialization of String?

/**
.Class.SparseString:
..cat:Index
..summary:A string storing only a fraction of the values of the original string..
..signature:SparseString<TValueString, TSpec>
..param.TValueString:The string containing the values.
..param.TSpec:The specialisation tag.
...default:void.
..include:seqan/String.h
*/
template <typename TValueString, typename TSpec = void>
struct SparseString
{
    typedef typename Fibre<SparseString, FibreValues>::Type         TFibreValues_;
    typedef typename Fibre<SparseString, FibreIndicators>::Type     TFibreIndicators_;
    typedef typename Size<SparseString>::Type                       TSize;

    TFibreValues_           values;
    TFibreIndicators_       indicators;
    TSize                   _length;

    SparseString() :
        _length(0)
    {};
};

// ==========================================================================
// Functions
// ==========================================================================

/*
.Function._assignValueInValueString
..param.container:
...type:Class.CompressedSA
*/
template <typename TFibreValues, typename TSpec, typename TPos, typename TValue>
inline void _assignValueInValueString(SparseString<TFibreValues, TSpec> & string, TPos pos, TValue value)
{
    getFibre(string, FibreValues())[pos] = value;
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

///.Function.clear.param.object.type:Class.SparseString
template <typename TFibreValues, typename TSpec>
inline void clear(SparseString<TFibreValues, TSpec> & string)
{
    string._length = 0;
    clear(getFibre(string, FibreValues()));
    clear(getFibre(string, FibreIndicators()));
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

///.Function.empty.param.object.type:Class.SparseString
template <typename TFibreValues, typename TSpec>
inline bool empty(SparseString<TFibreValues, TSpec> const & string)
{
//    return empty(getFibre(string, FibreIndicators()));
    return length(string) == 0;
}

template <typename TFibreValues, typename TSpec, typename TPos>
inline bool _isContained(SparseString<TFibreValues, TSpec> const & string, TPos const & pos)
{
    return getValue(getFibre(string, FibreIndicators()), pos);
}

// ----------------------------------------------------------------------------
// Function assignValue()
// ----------------------------------------------------------------------------

///.Function.assignValue.param.container.type:Class.SparseString
template <typename TFibreValues, typename TSpec, typename TPos, typename TValue>
inline void
assignValue(SparseString<TFibreValues, TSpec> & string, TPos pos, TValue value)
{
    if (!_isContained(string, pos))
        setValue(getFibre(string, FibreIndicators()), pos, false);

    getFibre(string, FibreValues())[getRank(getFibre(string, FibreIndicators()), pos) - 1] = value;
}

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------

///.Function.getValue.param.container.type:Class.SparseString
template <typename TFibreValues, typename TSpec, typename TPos>
inline typename GetValue<SparseString<TFibreValues, TSpec> >::Type
getValue(SparseString<TFibreValues, TSpec> & string, TPos pos)
{
    if (_isContained(string, pos))
        return getValue(getFibre(string, FibreValues()), getRank(getFibre(string, FibreIndicators()), pos) - 1);
    else 
        return DefaultValue<SparseString<TFibreValues, TSpec> >::VALUE;
}

template <typename TFibreValues, typename TSpec, typename TPos>
inline typename GetValue<SparseString<TFibreValues, TSpec> const>::Type
getValue(SparseString<TFibreValues, TSpec> const & string, TPos pos)
{
    if (_isContained(string, pos))
        return getValue(getFibre(string, FibreValues()), getRank(getFibre(string, FibreIndicators()), pos) - 1);
    else 
        return DefaultValue<SparseString<TFibreValues, TSpec> const>::VALUE;
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

///.Function.value.param.container.type:Class.SparseString
template <typename TFibreValues, typename TSpec, typename TPos>
inline typename Reference<SparseString<TFibreValues, TSpec> >::Type 
value(SparseString<TFibreValues, TSpec>&string, TPos pos)
{
    return getValue(string, pos);
}

template <typename TFibreValues, typename TSpec, typename TPos>
inline typename Reference<SparseString<TFibreValues, TSpec> >::Type
value(SparseString<TFibreValues, TSpec> const & string, TPos pos)
{
    return getValue(string, pos);
}

// ----------------------------------------------------------------------------
// Function getFibre()
// ----------------------------------------------------------------------------

///.Function.getFibre.param.container.type:Class.CompressedSA
template <typename TFibreValues, typename TSpec>
inline typename Fibre<SparseString<TFibreValues, TSpec>, FibreValues>::Type const &
getFibre(SparseString<TFibreValues, TSpec> const & sparseString, FibreValues)
{
    return sparseString.values;
}

template <typename TFibreValues, typename TSpec>
inline typename Fibre<SparseString<TFibreValues, TSpec>, FibreValues>::Type &
getFibre(SparseString<TFibreValues, TSpec>&sparseString, FibreValues)
{
    return sparseString.values;
}

template <typename TFibreValues, typename TSpec>
inline typename Fibre<SparseString<TFibreValues, TSpec>, FibreIndicators>::Type const &
getFibre(SparseString<TFibreValues, TSpec> const & sparseString, FibreIndicators)
{
    return sparseString.indicators;
}

template <typename TFibreValues, typename TSpec>
inline typename Fibre<SparseString<TFibreValues, TSpec>, FibreIndicators>::Type &
getFibre(SparseString<TFibreValues, TSpec>&sparseString, FibreIndicators)
{
    return sparseString.indicators;
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

///.Function.length.param.object.type:Class.SparseString
template <typename TFibreValues, typename TSpec>
inline typename Size<typename Fibre<SparseString<TFibreValues, TSpec>, FibreValues>::Type>::Type
length(SparseString<TFibreValues, TSpec> const & string)
{
    return string._length;
}

// ----------------------------------------------------------------------------
// Function resize()
// ----------------------------------------------------------------------------
// NOTE(esiragusa): This version of resize() was now working, therefore it was commented out.

///.Function.resize.param.object.type:Class.SparseString
//template <typename TFibreValues, typename TSpec, typename TSize, typename TValue, typename TExpand>
//inline typename Size<typename Fibre<SparseString<TFibreValues, TSpec>, FibreValues>::Type>::Type
//resize(SparseString<TFibreValues, TSpec> & string, TSize size, TValue value, Tag<TExpand> const tag)
//{
//    if (value != DefaultValue<SparseString<TFibreValues, TSpec> >::VALUE)
//    {
//        if (length(string) < size)
//            resize(getFibre(string, FibreValues()), length(getFibre(string, FibreValues())) + (size - length(string)), value, tag);
//        else
//            resize(getFibre(string, FibreValues()), getRank(getFibre(string, FibreIndicators()), size), value, tag);
//
//        resize(getFibre(string, FibreIndicators()), size, true, tag);
//    }
//
//    string._length = size;
//
//    return resize(getFibre(string, FibreIndicators()), size, false);
//}

template <typename TFibreValues, typename TSpec, typename TSize, typename TExpand>
inline typename Size<typename Fibre<SparseString<TFibreValues, TSpec>, FibreValues>::Type>::Type
resize(SparseString<TFibreValues, TSpec> & string, TSize size, Tag<TExpand> tag)
{
//    return resize(string, size, DefaultValue<SparseString<TFibreValues, TSpec> >::VALUE, tag);

    string._length = size;
    return resize(getFibre(string, FibreIndicators()), size, tag);
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

/**
.Function.open
..param.string:
...type:Class.SparseString
*/
template <typename TFibreValues, typename TSpec>
inline bool open(SparseString<TFibreValues, TSpec> & sparseString, const char * fileName, int openMode)
{
    String<char> name;

    // Length saved inside a .len file.
    name = fileName;
    append(name, ".len");
    if (!open(sparseString._length, toCString(name), openMode)) return false;

    // Values saved inside a .val file.
    name = fileName;
    append(name, ".val");
    if (!open(getFibre(sparseString, FibreValues()), toCString(name), openMode)) return false;

    // Indicators saved inside a .ind file.
    name = fileName;
    append(name, ".ind");
    if (!open(getFibre(sparseString, FibreIndicators()), toCString(name), openMode)) return false;

    return true;
}

template <typename TFibreValues, typename TSpec>
inline bool open(SparseString<TFibreValues, TSpec> & sparseString, const char * fileName)
{
    return open(sparseString, fileName, DefaultOpenMode<SparseString<TFibreValues, TSpec> >::VALUE);
}

// ----------------------------------------------------------------------------
// Function save()
// ----------------------------------------------------------------------------

/**
.Function.SparseString#save
..class:Class.SparseString
..summary:This functions saves a sparse string to disk.
..signature:open(string, fileName [, openMode])
..param.string:The string to be saved.
...type:Class.SparseString
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/
template <typename TFibreValues, typename TSpec>
inline bool save(SparseString<TFibreValues, TSpec> const & sparseString, const char * fileName)
{
    return save(sparseString, fileName, DefaultOpenMode<SparseString<TFibreValues, TSpec> >::VALUE);
}

template <typename TFibreValues, typename TSpec>
inline bool save(SparseString<TFibreValues, TSpec> const & sparseString, const char * fileName, int openMode)
{
    String<char> name;

    name = fileName;
    append(name, ".len");
    if (!save(length(sparseString), toCString(name), openMode)) return false;

    name = fileName;
    append(name, ".val");
    if (!save(getFibre(sparseString, FibreValues()), toCString(name), openMode)) return false;

    name = fileName;
    append(name, ".ind");
    if (!save(getFibre(sparseString, FibreIndicators()), toCString(name), openMode)) return false;

    return true;
}
// TODO(singer): setValue function

}
#endif // INDEX_FM_SPARSE_STRING_H_
