/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "LongList.H"
#include "Ostream.H"
#include "token.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, Foam::label Offset>
void Foam::LongList<T, Offset>::writeEntry(Ostream& os) const
{
    if
    (
        size() &&
        token::compound::isCompound
        (
            "List<" + word(pTraits<T>::typeName) + '>'
        )
    )
    {
        os  << word("List<" + word(pTraits<T>::typeName) + '>') << " ";
    }

    os << *this;
}

template<class T, Foam::label Offset>
void Foam::LongList<T, Offset>::writeEntry
(
    const word& keyword,
    Ostream& os
) const
{
    os.writeKeyword(keyword);
    writeEntry(os);
    os << token::END_STATEMENT << endl;
}

template<class T, Foam::label Offset>
Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::LongList<T, Offset>& DL
)
{
    if( (os.format() == IOstream::ASCII) || !contiguous<T>() )
    {
        if( DL.size() < 15 )
        {
            // Write size of list and start contents delimiter
            os << label(DL.size()) << token::BEGIN_LIST;

            // Write list contents
            forAll(DL, i)
            {
                if( i != 0 ) os << token::SPACE;
                os << DL[i];
            }

            // Write end of contents delimiter
            os << token::END_LIST;
        }
        else
        {
            // Write size of list and start contents delimiter
            os << nl << label(DL.size()) << nl << token::BEGIN_LIST;

            // Write list contents
            forAll(DL, i)
            {
                os << nl << DL[i];
            }

            // Write end of contents delimiter
            os << nl << token::END_LIST << nl;
        }
    }
    else
    {
        //- create a buffer list
        List<T> buf(DL.size());

        //- copy the data into the buffer
        forAll(buf, i)
            buf[i] = DL[i];

        //- write size
        os << nl << label(DL.size());

        //- write the buffer to stream
        os.write
        (
            reinterpret_cast<const char*>(buf.data()),
            label(buf.size() * sizeof(T))
        );

        os.flush();
    }

    // Check state of IOstream
    os.check("Ostream& operator<<(Ostream&, const LongList&)");

    return os;
}


template<class T, Foam::label Offset>
Foam::Istream& Foam::operator>>
(
    Foam::Istream& is,
    Foam::LongList<T, Offset>& DL
)
{
    // Anull list
    DL.setSize(0);

    is.fatalCheck("operator>>(Istream&, LongList<T, Offset>&)");

    token firstToken(is);

    is.fatalCheck
    (
        "operator>>(Istream&, LongList<T, Offset>&) : reading first token"
    );

    if( firstToken.isLabel() )
    {
        const sizeType size = firstToken.labelToken();

        // Set list length to that read
        DL.setSize(size);

        // Read list contents depending on data format
        if( (is.format() == IOstream::ASCII) || !contiguous<T>() )
        {
            // Read beginning of contents
            char listDelimiter = is.readBeginList("List");

            if( size == 0 )
            {
                if( listDelimiter != token::BEGIN_LIST )
                {
                    WarningIn
                    (
                        "template<class T, Foam::label Offset>"

                        "Foam::Istream& Foam::operator>>"
                        "("
                            "Foam::Istream& ,"
                            "Foam::LongList<T, Offset>& DL"
                        ")"
                    ) << "Missing ( after 0" << endl;

                    return is;
                }

                listDelimiter = is.readEndList("List");
                if( listDelimiter != token::END_LIST )
                {
                    WarningIn
                    (
                        "template<class T, Foam::label Offset>"

                        "Foam::Istream& Foam::operator>>"
                        "("
                            "Foam::Istream& ,"
                            "Foam::LongList<T, Offset>& DL"
                        ")"
                    ) << "Missing ) after 0(" << endl;
                }

                return is;
            }

            if( listDelimiter == token::BEGIN_LIST )
            {
                for(sizeType i=0;i<size;++i)
                {
                    is >> DL[i];

                    is.fatalCheck
                    (
                        "operator>>(Istream&, List<T>&) : reading entry"
                    );
                }
            }
            else
            {
                T element;
                is >> element;

                is.fatalCheck
                (
                    "operator>>(Istream&, List<T>&) : "
                    "reading the single entry"
                );

                for(sizeType i=0;i<size;++i)
                {
                    DL[i] = element;
                }
            }

            // Read end of contents
            is.readEndList("List");
        }
        else
        {
            List<T> buf(size);
            is.read(reinterpret_cast<char*>(buf.begin()), size * sizeof(T));

            forAll(buf, i)
                DL[i] = buf[i];

            buf.clear();

            is.fatalCheck
            (
                "operator>>(Istream&, LongList<T, Offset>&)"
                ": reading the binary block"
            );
        }
    }
    else
    {
        FatalErrorIn("operator>>(Istream&, LongList<T, Offset>&)")
            << "incorrect first token, expected <int>, found "
            << firstToken.info()
            << abort(FatalError);
    }

    return is;
}

template<class T, Foam::label Offset>
void Foam::LongList<T, Offset>::appendFromStream(Istream& is)
{
    is.fatalCheck("appendFromStream(Istream& is)");

    token firstToken(is);

    is.fatalCheck
    (
        "appendFromStream(Istream& is) : reading first token"
    );

    if( firstToken.isLabel() )
    {
        const sizeType size = firstToken.labelToken();

        if( size == 0 )
            return;

        sizeType origSize(this->size());

        // Set list length to that read
        setSize(origSize+size);

        // Read list contents depending on data format
        if( (is.format() == IOstream::ASCII) || !contiguous<T>() )
        {
            // Read beginning of contents
            char listDelimiter = is.readBeginList("List");

            if( listDelimiter == token::BEGIN_LIST )
            {
                for(sizeType i=0;i<size;++i)
                {
                    is >> this->operator[](origSize);
                    ++origSize;

                    is.fatalCheck
                    (
                        "appendFromStream(Istream& is) : reading entry"
                    );
                }
            }
            else
            {
                T element;
                is >> element;

                is.fatalCheck
                (
                    "appendFromStream(Istream& is) : "
                    "reading the single entry"
                );

                for(sizeType i=0;i<size;++i)
                {
                    this->operator[](origSize) = element;
                    ++origSize;
                }
            }

            // Read end of contents
            is.readEndList("List");
        }
        else
        {
            List<T> buf(size);
            is.read(reinterpret_cast<char*>(buf.begin()), size * sizeof(T));

            forAll(buf, i)
                this->operator[](origSize++) = buf[i];

            is.fatalCheck
            (
                "appendFromStream(Istream& is)"
                ": reading the binary block"
            );
        }
    }
    else
    {
        FatalErrorIn("appendFromStream(Istream& is)")
            << "incorrect first token, expected <int>, found "
            << firstToken.info()
            << abort(FatalError);
    }
}


// ************************************************************************* //
