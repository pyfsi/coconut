#include "DxfFile.H"

#include <iostream>

namespace Foam
{

// ============================================================================
// class DxfEntityArc
// ============================================================================

// ----------------------------------------------------------------------------
// Public member functions
// ----------------------------------------------------------------------------

vector DxfEntity::Ax(const vector& n) const
{
    if( mag(n.x()) < 1.0/64.0 && mag(n.y()) < 1.0/64.0 )
    {
        return (vector(0., 1.0, 0.0) ^ n);
    }

    return (vector(0.0, 0.0, 1.0) ^ n);
}

vector DxfEntity::Ay(const vector& n) const
{
    const vector Ax = this->Ax(n);

    return (n ^ Ax);
}

// ============================================================================
// class DxfEntityArc
// ============================================================================

// ----------------------------------------------------------------------------
// Public member functions
// ----------------------------------------------------------------------------

void DxfEntityArc::Print() const
{
    std::cout << "--- ARC ---" << std::endl
              << "  centre:     "
              << " " << centre_.x()
              << " " << centre_.y()
              << " " << centre_.z() << std::endl
              << "  radius:     "
              << " " << radius_ << std::endl
              << "  start angle:"
              << " " << startAngle_ << std::endl
              << "  end angle:  "
              << " " << endAngle_ << std::endl
              << "  normal:     "
              << " " << normal_.x()
              << " " << normal_.y()
              << " " << normal_.z() << std::endl
              << std::endl;
}


// ============================================================================
// class DxfEntityCircle
// ============================================================================

// ----------------------------------------------------------------------------
// Public member functions
// ----------------------------------------------------------------------------

void DxfEntityCircle::Print() const
{
    std::cout << "--- CIRCLE ---" << std::endl
              << "  centre:"
              << " " << centre_.x()
              << " " << centre_.y()
              << " " << centre_.z() << std::endl
              << "  radius:"
              << " " << radius_ << std::endl
              << "  normal:"
              << " " << normal_.x()
              << " " << normal_.y()
              << " " << normal_.z() << std::endl
              << std::endl;
}


// ============================================================================
// class DxfEntityEllipse
// ============================================================================

// ----------------------------------------------------------------------------
// Public member functions
// ----------------------------------------------------------------------------

void DxfEntityEllipse::Print() const
{
    std::cout << "--- ELLIPSE ---" << std::endl
              << "  centre:     "
              << " " << centre_.x()
              << " " << centre_.y()
              << " " << centre_.z() << std::endl
              << "  end major:  "
              << " " << endMajor_.x()
              << " " << endMajor_.y()
              << " " << endMajor_.z() << std::endl
              << "  axis ratio: "
              << " " << axisRatio_ << std::endl
              << "  start angle:"
              << " " << startAngle_ << std::endl
              << "  end angle:  "
              << " " << endAngle_ << std::endl
              << "  normal:     "
              << " " << normal_.x()
              << " " << normal_.y()
              << " " << normal_.z() << std::endl
              << std::endl;
}


// ============================================================================
// class DxfEntityInsert
// ============================================================================

// ----------------------------------------------------------------------------
// Public member functions
// ----------------------------------------------------------------------------

void DxfEntityInsert::Print() const
{
    std::cout << "--- INSERT ---" << std::endl
              << "  block name: " << blockName_ << std::endl
              << "  x insertion point: " << xInsertionPoint_ << std::endl
              << "  y insertion point: " << yInsertionPoint_ << std::endl
              << "  z insertion point: " << zInsertionPoint_ << std::endl
              << "  x scale factor: " << xScaleFactor_ << std::endl
              << "  y scale factor: " << yScaleFactor_ << std::endl
              << "  z scale factor: " << zScaleFactor_ << std::endl
              << "  rotation angle: " << rotationAngle_ << std::endl
              << "  column count: " << columnCount_ << std::endl
              << "  row count: " << rowCount_ << std::endl
              << "  column spacing: " << columnSpacing_ << std::endl
              << "  row spacing: " << rowSpacing_ << std::endl
              << "  x extrusion direction: " << xExtrusionDirection_ << std::endl
              << "  y extrusion direction: " << yExtrusionDirection_ << std::endl
              << "  z extrusion direction: " << zExtrusionDirection_ << std::endl
              << std::endl;
}


// ============================================================================
// class DxfEntityLine
// ============================================================================

// ----------------------------------------------------------------------------
// Public member functions
// ----------------------------------------------------------------------------

void DxfEntityLine::Print() const
{
    std::cout << "--- LINE ---" << std::endl
              << "  start: "
              << " " << start_.x()
              << " " << start_.y()
              << " " << start_.z() << std::endl
              << "  end:   "
              << " " << end_.x()
              << " " << end_.y()
              << " " << end_.z() << std::endl
              << "  normal:"
              << " " << normal_.x()
              << " " << normal_.y()
              << " " << normal_.z() << std::endl
              << std::endl;
}


// ============================================================================
// class DxfEntityPolyLine
// ============================================================================

// ----------------------------------------------------------------------------
// Public member functions
// ----------------------------------------------------------------------------

void DxfEntityPolyLine::Print() const
{
    std::cout << "--- LWPOLYLINE ---" << std::endl
              << std::endl;
}


// ============================================================================
// class DxfEntitySpline
// ============================================================================

// ----------------------------------------------------------------------------
// Public member functions
// ----------------------------------------------------------------------------

void DxfEntitySpline::Print() const
{
    std::cout << "--- SPLINE ---" << std::endl
              << std::endl;
}


// ============================================================================
// class DxfBlock
// ============================================================================

// ----------------------------------------------------------------------------
// Public member functions
// ----------------------------------------------------------------------------

void DxfBlock::Print() const
{
    std::cout << "--- BLOCK ---" << std::endl
              << "  " << name_ << std::endl;
    std::cout << "----------------------------------------------" << std::endl;

    for (auto it = entities_.begin(); it != entities_.end(); ++it)
    {
        (*it)->Print();
    }

    std::cout << "----------------------------------------------" << std::endl
              << std::endl;
}


// ============================================================================
// class DxfFile
// ============================================================================

// ----------------------------------------------------------------------------
// Public member functions
// ----------------------------------------------------------------------------

void DxfFile::Print() const
{
    std::cout << "--- UNITS ---" << std::endl
              << "  " << static_cast<int>(units_) << std::endl
              << std::endl;

    std::cout << "--- ANGLE UNITS ---" << std::endl
              << "  " << static_cast<int>(angleUnits_) << std::endl
              << std::endl;

    for (auto it = blocks_.begin(); it != blocks_.end(); ++it)
    {
        (*it)->Print();
    }

    for (auto it = entities_.begin(); it != entities_.end(); ++it)
    {
        (*it)->Print();
    }
}

} // End namespace Foam
