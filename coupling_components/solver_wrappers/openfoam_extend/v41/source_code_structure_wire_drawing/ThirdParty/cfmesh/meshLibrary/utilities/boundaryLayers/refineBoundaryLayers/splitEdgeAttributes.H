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

Class
    refineBoundaryLayers

Description
    Refine existing boundary layers

SourceFiles
    refineBoundaryLayers.C
    refineBoundaryLayersFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef splitEdgeAttributes_H
#define splitEdgeAttributes_H

#include "scalar.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class splitEdgeAttributes Declaration
\*---------------------------------------------------------------------------*/

class splitEdgeAttributes
{
    // Private data
        //- starting node label
        label start_;

        //- end node label
        label end_;

        //- number of splits
        label nSplits_;

        //- number of nodes at edge
        label nNodesAtEdge_;

        //- thickness of the first layer
        scalar firstLayerThickness_;

        //- thickness ratio
        scalar thicknessRatio_;

    public:

    // Constructors
        //- Null construct
        splitEdgeAttributes()
        :
            start_(-1),
            end_(),
            nSplits_(1),
            nNodesAtEdge_(0),
            firstLayerThickness_(-1.0),
            thicknessRatio_(1.0)
        {}

        //- Construct from components
        splitEdgeAttributes
        (
            const label start,
            const label end,
            const label nSplits,
            const label nNodesAtEdge,
            const scalar firstLayerThickness,
            const scalar thicknessRatio
        )
        :
            start_(start),
            end_(end),
            nSplits_(nSplits),
            nNodesAtEdge_(nNodesAtEdge),
            firstLayerThickness_(firstLayerThickness),
            thicknessRatio_(thicknessRatio)
        {}

    // Destructor
        ~splitEdgeAttributes()
        {}

    // Member functions
        //- return start node
        inline label start() const
        {
            return start_;
        }

        //- return end node
        inline label end() const
        {
            return end_;
        }

        //- return the number of splits
        inline label nSplits() const
        {
            return nSplits_;
        }

        //- return nodes at edge
        inline label nNodesAtEdge() const
        {
            return nNodesAtEdge_;
        }

        //- return first layer thickness
        inline scalar firstLayerThickness() const
        {
            return firstLayerThickness_;
        }

        //- return thickness ratio
        inline scalar thicknessRatio() const
        {
            return thicknessRatio_;
        }

    // Member operators

        inline void operator=(const splitEdgeAttributes& e)
        {
            start_ = e.start_;
            end_ = e.end_;
            nSplits_ = e.nSplits_;
            nNodesAtEdge_ = e.nNodesAtEdge_;
            firstLayerThickness_ = e.firstLayerThickness_;
            thicknessRatio_ = e.thicknessRatio_;
        }

        inline bool operator==(const splitEdgeAttributes& e) const
        {
            if( (start_ != e.start_) || (end_ != e.end_) )
                return false;

            return true;
        }

        inline bool operator!=(const splitEdgeAttributes& e) const
        {
            return !this->operator==(e);
        }

        inline bool operator<(const splitEdgeAttributes& e) const
        {
            if( start_ < e.start_ )
            {
                return true;
            }
            else if( (start_ == e.start_) && (end_ > e.end_) )
            {
                return true;
            }

            return false;
        }

    // Friend operators
        friend Ostream& operator<<
        (
            Ostream& os,
            const splitEdgeAttributes& e
        )
        {
            os << token::BEGIN_LIST;
            os << e.start_ << token::SPACE;
            os << e.end_ << token::SPACE << token::BEGIN_LIST;
            os << e.nSplits_ << token::SPACE;
            os << e.nNodesAtEdge_ << token::SPACE;
            os << e.firstLayerThickness_ << token::SPACE;
            os << e.thicknessRatio_ << token::END_LIST;
            os << token::END_LIST;

            // Check state of Ostream
            os.check
            (
                "operator<<(Ostream&, const splitEdgeAttributes&"
            );

            return os;
        }

        friend Istream& operator>>(Istream& is, splitEdgeAttributes& e)
        {
            // Read beginning of splitEdgeAttributes
            is.readBegin("splitEdgeAttributes");

            is >> e.start_;
            is >> e.end_;
            is >> e.nSplits_;
            is >> e.nNodesAtEdge_;
            is >> e.firstLayerThickness_;
            is >> e.thicknessRatio_;

            // Read end of labelledPair
            is.readEnd("splitEdgeAttributes");

            // Check state of Istream
            is.check("operator>>(Istream&, splitEdgeAttributes");

            return is;
        }
};

//- Specify data associated with splitEdgeAttributes
//- type is contiguous
template<>
inline bool contiguous<splitEdgeAttributes>()
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
