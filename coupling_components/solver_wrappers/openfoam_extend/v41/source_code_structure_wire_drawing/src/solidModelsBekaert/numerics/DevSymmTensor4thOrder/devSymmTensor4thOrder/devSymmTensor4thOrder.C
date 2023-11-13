/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "devSymmTensor4thOrder.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{



    
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const devSymmTensor4thOrder::typeName = "devSymmTensor4thOrder";

template<>
const char* devSymmTensor4thOrder::componentNames[] =
{
     	  "xxxy",
	  "yzyz",
	  "xxyy",
	  "yyyz",
	  "xzyy",
	  "xyxy",
	  "xyyy",
	  "xxxz",
	  "xxyz",
	  "xzyz",
	  "yyyy",
	  "xxxx",
	  "xyyz",
	  "xzxz",
	  "xyxz"
};

template<>
const devSymmTensor4thOrder devSymmTensor4thOrder::zero
(
    0.0,
	  0.0,
	  0.0,
	  0.0,
	  0.0,
	  0.0,
	  0.0,
	  0.0,
	  0.0,
	  0.0,
	  0.0,
	  0.0,
	  0.0,
	  0.0,
	  0.0
);

// template<>
// const devSymmTensor4thOrder devSymmTensor4thOrder::one
// (
//      	  0.0,
// 	  3.0/4.0*(2.0/3.0),
// 	  -0.5*(2.0/3.0),
// 	  0.0,
// 	  0.0,
// 	  3.0/4.0*(2.0/3.0),
// 	  0.0,
// 	  0.0,
// 	  0.0,
// 	  0.0,
// 	  1.0*(2.0/3.0),
// 	  1.0*(2.0/3.0),
// 	  0.0,
// 	  3.0/4.0*(2.0/3.0),
// 	  0.0
// );
template<>
const devSymmTensor4thOrder devSymmTensor4thOrder::one
(
     	  0.0,
	  3.0/4.0,
	  -0.5,
	  0.0,
	  0.0,
	  3.0/4.0,
	  0.0,
	  0.0,
	  0.0,
	  0.0,
	  1.0,
	  1.0,
	  0.0,
	  3.0/4.0,
	  0.0
);

template<>
const devSymmTensor4thOrder devSymmTensor4thOrder::max
(
     	  VGREAT,
	  VGREAT,
	  VGREAT,
	  VGREAT,
	  VGREAT,
	  VGREAT,
	  VGREAT,
	  VGREAT,
	  VGREAT,
	  VGREAT,
	  VGREAT,
	  VGREAT,
	  VGREAT,
	  VGREAT,
	  VGREAT
);

template<>
const devSymmTensor4thOrder devSymmTensor4thOrder::min
(
     	  -VGREAT,
	  -VGREAT,
	  -VGREAT,
	  -VGREAT,
	  -VGREAT,
	  -VGREAT,
	  -VGREAT,
	  -VGREAT,
	  -VGREAT,
	  -VGREAT,
	  -VGREAT,
	  -VGREAT,
	  -VGREAT,
	  -VGREAT,
	  -VGREAT
);




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
