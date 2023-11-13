#ifndef DxfFileToTriSurfConverter_I_H
#define DxfFileToTriSurfConverter_I_H

namespace Foam
{

// ============================================================================
// class DxfFileToTriSurfConverter
// ============================================================================

// ----------------------------------------------------------------------------
// Public member functions
// ----------------------------------------------------------------------------

template <class InputIt>
void DxfFileToTriSurfConverter::AddEntityToTriSurf
(
    InputIt entity,
    triSurf& surf
)
{
    if (typeid(**entity) == typeid(DxfEntityArc))
    {
        AddArcToTriSurf
        (
            dynamic_cast<const DxfEntityArc&>(**entity),
            surf
        );
    }
    else if (typeid(**entity) == typeid(DxfEntityCircle))
    {
        AddCircleToTriSurf
        (
            dynamic_cast<const DxfEntityCircle&>(**entity),
            surf
        );
    }
    else if (typeid(**entity) == typeid(DxfEntityEllipse))
    {
        AddEllipseToTriSurf
        (
            dynamic_cast<const DxfEntityEllipse&>(**entity),
            surf
        );
    }
    else if (typeid(**entity) == typeid(DxfEntityInsert))
    {
        AddInsertToTriSurf
        (
            dynamic_cast<const DxfEntityInsert&>(**entity),
            surf
        );
    }
    else if (typeid(**entity) == typeid(DxfEntityLine))
    {
        AddLineToTriSurf
        (
            dynamic_cast<const DxfEntityLine&>(**entity),
            surf
        );
    }
    else if (typeid(**entity) == typeid(DxfEntityPolyLine))
    {
        AddPolyLineToTriSurf
        (
            dynamic_cast<const DxfEntityPolyLine&>(**entity),
            surf
        );
    }
    else if (typeid(**entity) == typeid(DxfEntitySpline))
    {
        AddSplineToTriSurf
        (
            dynamic_cast<const DxfEntitySpline&>(**entity),
            surf
        );
    }
    else
    {
        throw DxfFileToTriSurfException("Unknown ENTITY");
    }
}

} // End namespace Foam

#endif // DxfFileToTriSurfConverter_I_H
