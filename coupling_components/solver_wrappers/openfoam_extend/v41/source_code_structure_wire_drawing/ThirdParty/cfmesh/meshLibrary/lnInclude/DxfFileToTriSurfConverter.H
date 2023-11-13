#ifndef DxfFileToTriSurfConverter_H
#define DxfFileToTriSurfConverter_H

#include "DxfFile.H"
#include "scalar.H"

#include <stdexcept>

namespace Foam
{

// ============================================================================
// Predeclarations
// ============================================================================

class triSurf;


// ============================================================================
// class DxfFileToTriSurfException
// ============================================================================

class DxfFileToTriSurfException : public std::runtime_error
{
    using std::runtime_error::runtime_error;
};


// ============================================================================
// class DxfFileToTriSurfConverter
// ============================================================================

class DxfFileToTriSurfConverter
{
public:

    // ------------------------------------------------------------------------
    // Public member functions
    // ------------------------------------------------------------------------

    DxfFileToTriSurfConverter
    (
        const DxfFile& file,
        const scalar&  angleTol = 1.0,
        const scalar&  discretisationLength = VGREAT
    );

    void Convert
    (
        triSurf& surf
    );

private:

    // ------------------------------------------------------------------------
    // Private member functions
    // ------------------------------------------------------------------------

    template <class InputIt>
    void AddEntityToTriSurf
    (
        InputIt entity,
        triSurf& surf
    );

    void AddArcToTriSurf
    (
        const DxfEntityArc& arc,
        triSurf& surf
    );

    void AddCircleToTriSurf
    (
        const DxfEntityCircle& circle,
        triSurf& surf
    );

    void AddEllipseToTriSurf
    (
        const DxfEntityEllipse& ellipse,
        triSurf& surf
    );

    void AddInsertToTriSurf
    (
        const DxfEntityInsert& insert,
        triSurf& surf
    );

    void AddLineToTriSurf
    (
        const DxfEntityLine& line,
        triSurf& surf
    );

    void AddPolyLineToTriSurf
    (
        const DxfEntityPolyLine& polyLine,
        triSurf& surf
    );

    void AddSplineToTriSurf
    (
        const DxfEntitySpline& spline,
        triSurf& surf
    );

    void ScaleToMeters(triSurf& surf);

    // ------------------------------------------------------------------------
    // Private data members
    // ------------------------------------------------------------------------

    const DxfFile& file_;

    scalar angleTol_;
    scalar discretisationLength_;

};

} // End namespace Foam

#include "DxfFileToTriSurfConverter.I.H"

#endif // DxfFileToTriSurfConverter_H
