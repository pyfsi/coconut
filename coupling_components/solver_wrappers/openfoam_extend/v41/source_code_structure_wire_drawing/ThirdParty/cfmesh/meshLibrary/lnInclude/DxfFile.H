#ifndef DxfFile_H
#define DxfFile_H

#include "point.H"
#include "vector.H"

#include <list>
#include <memory>
#include <string>
#include <vector>

namespace Foam
{

// ============================================================================
// class DxfEntity
// ============================================================================

class DxfEntity
{
public:

    // ------------------------------------------------------------------------
    // Public member functions
    // ------------------------------------------------------------------------

    virtual ~DxfEntity()
    { }

    virtual void Print() const
    { }

    vector Ax(const vector& n) const;
    vector Ay(const vector& n) const;
};


// ============================================================================
// class DxfEntityArc
// ============================================================================

class DxfEntityArc : public DxfEntity
{
public:

    // ------------------------------------------------------------------------
    // Public member functions
    // ------------------------------------------------------------------------

    DxfEntityArc
    (
        const point&  centre,
        const double& radius,
        const double& startAngle,
        const double& endAngle,
        const vector& normal
    )
        : centre_(centre)
        , radius_(radius)
        , startAngle_(startAngle)
        , endAngle_(endAngle)
        , normal_(normal)
    { }

    const point&  GetCentre() const     { return centre_; }
    const double& GetRadius() const     { return radius_; }
    const double& GetStartAngle() const { return startAngle_; }
    const double& GetEndAngle() const   { return endAngle_; }
    const vector& GetNormal() const     { return normal_; }

    void Print() const override;

private:

    // ------------------------------------------------------------------------
    // Private data members
    // ------------------------------------------------------------------------

    point  centre_;
    double radius_;
    double startAngle_;
    double endAngle_;
    vector normal_;

};


// ============================================================================
// class DxfEntityCircle
// ============================================================================

class DxfEntityCircle : public DxfEntity
{
public:

    // ------------------------------------------------------------------------
    // Public member functions
    // ------------------------------------------------------------------------

    DxfEntityCircle
    (
        const point&  centre,
        const double& radius,
        const vector& normal
    )
        : centre_(centre)
        , radius_(radius)
        , normal_(normal)
    { }

    const point&  GetCentre() const { return centre_; }
    const double& GetRadius() const { return radius_; }
    const vector& GetNormal() const { return normal_; }

    void Print() const override;

private:

    // ------------------------------------------------------------------------
    // Private data members
    // ------------------------------------------------------------------------

    point  centre_;
    double radius_;
    vector normal_;

};


// ============================================================================
// class DxfEntityEllipse
// ============================================================================

class DxfEntityEllipse : public DxfEntity
{
public:

    // ------------------------------------------------------------------------
    // Public member functions
    // ------------------------------------------------------------------------

    DxfEntityEllipse
    (
        const point&  centre,
        const point&  endMajor,
        const vector& normal,
        const double& axisRatio,
        const double& startAngle,
        const double& endAngle
    )
        : centre_(centre)
        , endMajor_(endMajor)
        , normal_(normal)
        , axisRatio_(axisRatio)
        , startAngle_(startAngle)
        , endAngle_(endAngle)
    { }

    const point&  GetCentre() const     { return centre_; }
    const point&  GetEndMajor() const   { return endMajor_; }
    const vector& GetNormal() const     { return normal_; }
    const double& GetAxisRatio() const  { return axisRatio_; }
    const double& GetStartAngle() const { return startAngle_; }
    const double& GetEndAngle() const   { return endAngle_; }

    void Print() const override;

public:

    // ------------------------------------------------------------------------
    // Private data members
    // ------------------------------------------------------------------------

    point  centre_;
    point  endMajor_;
    vector normal_;
    double axisRatio_;
    double startAngle_;
    double endAngle_;

};


// ============================================================================
// class DxfEntityInsert
// ============================================================================

class DxfEntityInsert : public DxfEntity
{
public:

    // ------------------------------------------------------------------------
    // Public member functions
    // ------------------------------------------------------------------------

    DxfEntityInsert
    (
        const std::string blockName,

        const double xInsertionPoint,
        const double yInsertionPoint,
        const double zInsertionPoint,

        const double xScaleFactor,
        const double yScaleFactor,
        const double zScaleFactor,

        const double rotationAngle,

        const int columnCount,
        const int rowCount,

        const double columnSpacing,
        const double rowSpacing,

        const double xExtrusionDirection,
        const double yExtrusionDirection,
        const double zExtrusionDirection
    )
        : blockName_(blockName)

        , xInsertionPoint_(xInsertionPoint)
        , yInsertionPoint_(yInsertionPoint)
        , zInsertionPoint_(zInsertionPoint)

        , xScaleFactor_(xScaleFactor)
        , yScaleFactor_(yScaleFactor)
        , zScaleFactor_(zScaleFactor)

        , rotationAngle_(rotationAngle)

        , columnCount_(columnCount)
        , rowCount_(rowCount)

        , columnSpacing_(columnSpacing)
        , rowSpacing_(rowSpacing)

        , xExtrusionDirection_(xExtrusionDirection)
        , yExtrusionDirection_(yExtrusionDirection)
        , zExtrusionDirection_(zExtrusionDirection)
    { }

    const std::string& GetBlockName() const { return blockName_; }

    void Print() const override;

private:

    // ------------------------------------------------------------------------
    // Private data members
    // ------------------------------------------------------------------------

    std::string blockName_;

    double xInsertionPoint_;
    double yInsertionPoint_;
    double zInsertionPoint_;

    double xScaleFactor_;
    double yScaleFactor_;
    double zScaleFactor_;

    double rotationAngle_;

    int columnCount_;
    int rowCount_;

    double columnSpacing_;
    double rowSpacing_;

    double xExtrusionDirection_;
    double yExtrusionDirection_;
    double zExtrusionDirection_;

};


// ============================================================================
// class DxfEntityLine
// ============================================================================

class DxfEntityLine : public DxfEntity
{
public:

    // ------------------------------------------------------------------------
    // Public member functions
    // ------------------------------------------------------------------------

    DxfEntityLine
    (
        const point&  start,
        const point&  end,
        const vector& normal
    )
        : start_(start)
        , end_(end)
        , normal_(normal)
    { }

    const point&  GetStart() const  { return start_; }
    const point&  GetEnd() const    { return end_; }
    const vector& GetNormal() const { return normal_; }

    void Print() const override;

private:

    // ------------------------------------------------------------------------
    // Private data members
    // ------------------------------------------------------------------------

    point  start_;
    point  end_;
    vector normal_;

};


// ============================================================================
// class DxfEntityPolyLine
// ============================================================================

class DxfEntityPolyLine : public DxfEntity
{
public:

    // ------------------------------------------------------------------------
    // Public member functions
    // ------------------------------------------------------------------------

    DxfEntityPolyLine
    (
        const int& numOfVertices,
        const int& flag,
        const std::vector<double>& verticesX,
        const std::vector<double>& verticesY,
        const vector& normal,
        const scalar elevation
    )
        : numOfVertices_(numOfVertices)
        , flag_(flag)
        , verticesX_(verticesX)
        , verticesY_(verticesY)
        , normal_(normal)
        , elevation_(elevation)
    { }

    const int& GetNumOfVertices() const { return numOfVertices_; }
    const int& GetFlag() const          { return flag_; }

    const std::vector<double>& GetVerticesX() const { return verticesX_; }
    const std::vector<double>& GetVerticesY() const { return verticesY_; }

    const vector& GetNormal() const { return normal_; }

    const scalar& GetElevation() const { return elevation_; }

    void Print() const override;

private:

    // ------------------------------------------------------------------------
    // Private data members
    // ------------------------------------------------------------------------

    int numOfVertices_;
    int flag_;
    std::vector<double> verticesX_;
    std::vector<double> verticesY_;
    vector normal_;
    scalar elevation_;
};


// ============================================================================
// class DxfEntitySpline
// ============================================================================

class DxfEntitySpline : public DxfEntity
{
public:

    // ------------------------------------------------------------------------
    // Public member functions
    // ------------------------------------------------------------------------

    DxfEntitySpline
    (
        const vector& normal,
        const int&    flag,
        const int&    degree,
        const int&    numOfKnots,
        const int&    numOfControlPoints,
        const int&    numOfFitPoints,
        const double& knotTol,
        const double& controlPointTol,
        const double& fitTol,
        const vector& startTangent,
        const vector& endTangent,
        const std::list<double>& knotValues,
        const std::list<double>& weights,
        const std::list<double>& controlPointsX,
        const std::list<double>& controlPointsY,
        const std::list<double>& controlPointsZ,
        const std::list<double>& fitPointsX,
        const std::list<double>& fitPointsY,
        const std::list<double>& fitPointsZ
    )
        : normal_(normal)
        , flag_(flag)
        , degree_(degree)
        , numOfKnots_(numOfKnots)
        , numOfControlPoints_(numOfControlPoints)
        , numOfFitPoints_(numOfFitPoints)
        , knotTol_(knotTol)
        , controlPointTol_(controlPointTol)
        , fitTol_(fitTol)
        , startTangent_(startTangent)
        , endTangent_(endTangent)
        , knotValues_(knotValues)
        , weights_(weights)
        , controlPointsX_(controlPointsX)
        , controlPointsY_(controlPointsY)
        , controlPointsZ_(controlPointsZ)
        , fitPointsX_(fitPointsX)
        , fitPointsY_(fitPointsY)
        , fitPointsZ_(fitPointsZ)
    { }

    int GetFlag() const
    {
        return flag_;
    }

    int GetNumOfControlPoints() const
    {
        return numOfControlPoints_;
    }

    const vector& GetStartTangent() const
    {
        return startTangent_;
    }

    const vector& GetEndTangent() const
    {
        return endTangent_;
    }

    const std::list<double>& GetControlPointsX() const
    {
        return controlPointsX_;
    }

    const std::list<double>& GetControlPointsY() const
    {
        return controlPointsY_;
    }

    const std::list<double>& GetControlPointsZ() const
    {
        return controlPointsZ_;
    }

    void Print() const override;

private:

    // ------------------------------------------------------------------------
    // Private data members
    // ------------------------------------------------------------------------

    vector normal_;
    int    flag_;
    int    degree_;
    int    numOfKnots_;
    int    numOfControlPoints_;
    int    numOfFitPoints_;
    double knotTol_ = 0.0000001;
    double controlPointTol_ = 0.0000001;
    double fitTol_ = 0.0000000001;
    vector startTangent_;
    vector endTangent_;

    std::list<double> knotValues_;
    std::list<double> weights_;

    std::list<double> controlPointsX_;
    std::list<double> controlPointsY_;
    std::list<double> controlPointsZ_;

    std::list<double> fitPointsX_;
    std::list<double> fitPointsY_;
    std::list<double> fitPointsZ_;

};


// ============================================================================
// class DxfBlock
// ============================================================================

class DxfBlock
{
public:

    // ------------------------------------------------------------------------
    // Public member types
    // ------------------------------------------------------------------------

    using EntityList = std::list<std::shared_ptr<DxfEntity>>;

    // ------------------------------------------------------------------------
    // Public member functions
    // ------------------------------------------------------------------------

    DxfBlock(const std::string& name)
        : name_(name)
    { }

    const std::string& GetName() const
    {
        return name_;
    }

    void AddEntity(std::shared_ptr<DxfEntity> entity)
    {
        entities_.push_back(entity);
    }

    const EntityList& GetEntities() const
    {
        return entities_;
    }

    void Print() const;

private:

    // ------------------------------------------------------------------------
    // Private data members
    // ------------------------------------------------------------------------

    std::string name_;
    EntityList entities_;

};


// ============================================================================
// class DxfFile
// ============================================================================

class DxfFile
{
public:

    // ------------------------------------------------------------------------
    // Public member types
    // ------------------------------------------------------------------------

    enum class Units
    {
        Unitless = 0,
        Inches = 1,
        Feet = 2,
        Miles = 3,
        Millimeters = 4,
        Centimeters = 5,
        Meters = 6,
        Kilometers = 7,
        Microinches = 8,
        Mils = 9,
        Yards = 10,
        Angstroms = 11,
        Nanometers = 12,
        Microns = 13,
        Decimeters = 14,
        Decameters = 15,
        Hectometers = 16,
        Gigameters = 17,
        AstronomicalUnits = 18,
        LightYears = 19,
        Parsecs = 20
    };

    enum class AngleUnits
    {
        DecimalDegree = 0,
        DegMinSec = 1,
        Gradians = 2,
        Radians = 3,
        SurveyorsUnits = 4
    };

    using BlockList  = std::list<std::shared_ptr<DxfBlock>>;
    using EntityList = std::list<std::shared_ptr<DxfEntity>>;

    // ------------------------------------------------------------------------
    // Public member functions
    // ------------------------------------------------------------------------

    //
    // Units
    //

    Units GetUnits() const
    {
        return units_;
    }

    void SetUnits(Units units)
    {
        units_ = units;
    }

    //
    // Angle units
    //

    AngleUnits GetAngleUnits() const
    {
        return angleUnits_;
    }

    void SetAngleFormat(AngleUnits angleUnits)
    {
        angleUnits_ = angleUnits;
    }

    //
    // Blocks
    //

    void AddBlock(std::shared_ptr<DxfBlock> block)
    {
        blocks_.push_back(block);
    }

    const BlockList& GetBlocks() const
    {
        return blocks_;
    }

    //
    // Entities
    //

    void AddEntity(std::shared_ptr<DxfEntity> entity)
    {
        entities_.push_back(entity);
    }

    const EntityList& GetEntities() const
    {
        return entities_;
    }

    void Print() const;

private:

    // ------------------------------------------------------------------------
    // Private data members
    // ------------------------------------------------------------------------

    Units      units_      = Units::Unitless;
    AngleUnits angleUnits_ = AngleUnits::DecimalDegree;

    BlockList  blocks_;
    EntityList entities_;

};

} // End namespace Foam

#endif // DxfFile_H
