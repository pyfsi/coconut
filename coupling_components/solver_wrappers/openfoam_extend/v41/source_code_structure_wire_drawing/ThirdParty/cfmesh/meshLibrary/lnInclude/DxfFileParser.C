#include "DxfFileParser.H"

#include <iostream>
#include <string>
#include <sstream>

namespace Foam
{

// ============================================================================
// Non-member functions
// ============================================================================

std::istream& operator>>(std::istream& is, DxfRecord& record)
{
    std::istringstream iss;
    std::string line;

    // ------------------
    // Reading group code
    // ------------------

    std::getline(is, line);

    if (is.eof())
    {
        throw DxfFileParserException("Encountered EOF.");
    }
    else if (is.fail())
    {
        throw DxfFileParserException("Error while reading line.");
    }

    // Removing trailing <CR> and <LF> characters
    while (!line.empty() && (*line.rbegin() == '\r' || *line.rbegin() == '\n'))
    {
        line.erase(line.size() - 1, 1);
    }

    // Removing leading spaces and tabs
    while (!line.empty() && (*line.begin() == ' ' || *line.begin() == '\t'))
    {
        line.erase(0, 1);
    }

    // Checking if line is empty
    if (line.empty())
    {
        throw DxfFileParserException("Empty line. Expected group code.");
    }

    // Reading group code
    int code;
    iss.str(line);
    iss >> code;

    if (iss.fail())
    {
        throw DxfFileParserException("Error while reading group code");
    }

    record.SetCode(code);

    // -------------
    // Reading value
    // -------------

    std::getline(is, line);

    if (is.fail())
    {
        throw DxfFileParserException("Error while reading line.");
    }

    // Removing trailing <CR> and <LF> characters
    while (!line.empty() && (*line.rbegin() == '\r' || *line.rbegin() == '\n'))
    {
        line.erase(line.size() - 1, 1);
    }

    // Removing leading spaces and tabs
    while (!line.empty() && (*line.begin() == ' ' || *line.begin() == '\t'))
    {
        line.erase(0, 1);
    }

    record.SetValue(line);

    return is;
}


// ============================================================================
// class DxfFileParser
// ============================================================================

// ----------------------------------------------------------------------------
// Public member functions
// ----------------------------------------------------------------------------

void DxfFileParser::Parse(std::istream& is, DxfFile& file)
{
    IsItType isit(is);

    while (isit->GetCode() != 0 || isit->GetValue() != "EOF")
    {
        if (isit->GetCode() == 0 && isit->GetValue() == "SECTION")
        {
            ++isit;

            if (isit->GetCode() == 2 && isit->GetValue() == "HEADER")
            {
                ParseSectionHeader(isit, file);
            }
            else if (isit->GetCode() == 2 && isit->GetValue() == "BLOCKS")
            {
                ParseSectionBlocks(isit, file);
            }
            else if (isit->GetCode() == 2 && isit->GetValue() == "ENTITIES")
            {
                ParseSectionEntities(isit, file);
            }
        }

        ++isit;
    }
}


// ----------------------------------------------------------------------------
// Private member functions
// ----------------------------------------------------------------------------

void DxfFileParser::ParseSectionHeader
(
    DxfFileParser::IsItType& isit,
    DxfFile& file
)
{
    // Before incrementing isit points to (2, HEADER) record
    ++isit;

    while (isit->GetCode() != 0 || isit->GetValue() != "ENDSEC")
    {
        if (isit->GetCode() == 9 && isit->GetValue() == "$INSUNITS")
        {
            file.SetUnits(ParseUnits(isit));
        }
        else if (isit->GetCode() == 9 && isit->GetValue() == "$DIMAUNIT")
        {
            file.SetAngleFormat(ParseAngleUnits(isit));
        }

        ++isit;
    }
}


DxfFile::Units DxfFileParser::ParseUnits
(
    DxfFileParser::IsItType& isit
)
{
    // Before incrementing isit points to (9, $INSUNITS) record
    ++isit;

    if (isit->GetCode() != 70)
    {
        throw DxfFileParserException
        (
            "Unexpected group code while parsing $INSUNITS value."
        );
    }

    int units;
    std::istringstream iss(isit->GetValue());
    iss >> units;
    if (iss.fail())
    {
        throw DxfFileParserException
        (
            "Error while reading $INSUNITS value."
        );
    }

    return static_cast<DxfFile::Units>(units);
}


DxfFile::AngleUnits DxfFileParser::ParseAngleUnits
(
    DxfFileParser::IsItType& isit
)
{
    // Before incrementing isit points to (9, $DIMAUNIT) record
    ++isit;

    if (isit->GetCode() != 70)
    {
        throw DxfFileParserException
        (
            "Unexpected group code while parsing $DIMAUNIT value."
        );
    }

    int angleUnits;
    std::istringstream iss(isit->GetValue());
    iss >> angleUnits;
    if (iss.fail())
    {
        throw DxfFileParserException
        (
            "Error while reading $DIMAUNIT value."
        );
    }

    return static_cast<DxfFile::AngleUnits>(angleUnits);
}


void DxfFileParser::ParseSectionBlocks
(
    DxfFileParser::IsItType& isit,
    DxfFile& file
)
{
    // Before incrementing isit points to (2, BLOCKS) record
    ++isit;

    while (isit->GetCode() != 0 || isit->GetValue() != "ENDSEC")
    {
        if (isit->GetCode() == 0 && isit->GetValue() == "BLOCK")
        {
            file.AddBlock(ParseBlock(isit));
        }

        ++isit;
    }
}


std::shared_ptr<DxfBlock> DxfFileParser::ParseBlock
(
    DxfFileParser::IsItType& isit
)
{
    std::string blockName;
    bool        blockNameRead = false;

    // Before incrementing isit points to (0, BLOCK) record
    ++isit;

    // Reading block entity definition
    while (isit->GetCode() != 0)
    {
        if (isit->GetCode() == 2)
        {
            // Reading (2, <block name>) record
            std::istringstream iss(isit->GetValue());
            iss >> blockName;
            if (iss.fail())
            {
                throw DxfFileParserException
                (
                    "Error while reading <block name> value."
                );
            }
            blockNameRead = true;
        }

        ++isit;
    }

    if (!blockNameRead)
    {
        throw DxfFileParserException
        (
            "Cannot find record (2, <block name>)."
        );
    }

    auto block = std::make_shared<DxfBlock>(blockName);

    // Now isit points to some record with group code 0. Continue parsing
    // until ENDBLK is encountered
    while (isit->GetCode() != 0 || isit->GetValue() != "ENDBLK")
    {
        if (isit->GetCode() == 0 && isit->GetValue() == "ARC")
        {
            block->AddEntity(ParseArc(isit));
        }
        else if (isit->GetCode() == 0 && isit->GetValue() == "CIRCLE")
        {
            block->AddEntity(ParseCircle(isit));
        }
        else if (isit->GetCode() == 0 && isit->GetValue() == "ELLIPSE")
        {
            block->AddEntity(ParseEllipse(isit));
        }
        else if (isit->GetCode() == 0 && isit->GetValue() == "LINE")
        {
            block->AddEntity(ParseLine(isit));
        }
        else if (isit->GetCode() == 0 && isit->GetValue() == "LWPOLYLINE")
        {
            block->AddEntity(ParseLwPolyLine(isit));
        }
        else if (isit->GetCode() == 0 && isit->GetValue() == "SPLINE")
        {
            block->AddEntity(ParseSpline(isit));
        }
        else
        {
            ++isit;
        }
    }

    return block;
}


void DxfFileParser::ParseSectionEntities
(
    DxfFileParser::IsItType& isit,
    DxfFile& file
)
{
    ++isit; // Now isit points to first <entity type> record

    while (isit->GetCode() != 0 || isit->GetValue() != "ENDSEC")
    {
        if (isit->GetCode() == 0 && isit->GetValue() == "ARC")
        {
            file.AddEntity(ParseArc(isit));
        }
        else if (isit->GetCode() == 0 && isit->GetValue() == "CIRCLE")
        {
            file.AddEntity(ParseCircle(isit));
        }
        else if (isit->GetCode() == 0 && isit->GetValue() == "ELLIPSE")
        {
            file.AddEntity(ParseEllipse(isit));
        }
        else if (isit->GetCode() == 0 && isit->GetValue() == "INSERT")
        {
            file.AddEntity(ParseInsert(isit));
        }
        else if (isit->GetCode() == 0 && isit->GetValue() == "LINE")
        {
            file.AddEntity(ParseLine(isit));
        }
        else if (isit->GetCode() == 0 && isit->GetValue() == "LWPOLYLINE")
        {
            file.AddEntity(ParseLwPolyLine(isit));
        }
        else if (isit->GetCode() == 0 && isit->GetValue() == "SPLINE")
        {
            file.AddEntity(ParseSpline(isit));
        }
        else
        {
            ++isit;
        }
    }
}


std::shared_ptr<DxfEntity> DxfFileParser::ParseArc
(
    DxfFileParser::IsItType& isit
)
{
    double centreX = 0;
    double centreY = 0;
    double centreZ = 0;
    double radius  = 0;
    double startAngle = 0;
    double endAngle   = 0;
    double normalX = 0;
    double normalY = 0;
    double normalZ = 1;

    ++isit;
    while (isit->GetCode() != 0)
    {
        double value;

        switch (isit->GetCode())
        {
            case 10:
            case 20:
            case 30:
            case 40:
            case 50:
            case 51:
            case 210:
            case 220:
            case 230:
            {
                std::istringstream iss(isit->GetValue());
                iss >> value;
                if (iss.fail())
                {
                    throw DxfFileParserException("Error while reading value.");
                }
                break;
            }

            default:
                break;
        }

        switch (isit->GetCode())
        {
            case 10:  centreX = value;    break;
            case 20:  centreY = value;    break;
            case 30:  centreZ = value;    break;
            case 40:  radius = value;     break;
            case 50:  startAngle = value; break;
            case 51:  endAngle = value;   break;
            case 210: normalX = value;    break;
            case 220: normalY = value;    break;
            case 230: normalZ = value;    break;
            default: break;
        }

        ++isit;
    }

    return std::make_shared<DxfEntityArc>
    (
        point(centreX, centreY, centreZ),
        radius,
        startAngle,
        endAngle,
        vector(normalX, normalY, normalZ)
    );
}


std::shared_ptr<DxfEntity> DxfFileParser::ParseCircle
(
    DxfFileParser::IsItType& isit
)
{
    double centreX = 0;
    double centreY = 0;
    double centreZ = 0;
    double radius  = 0;
    double normalX = 0;
    double normalY = 0;
    double normalZ = 1;

    ++isit;
    while (isit->GetCode() != 0)
    {
        double value;

        switch (isit->GetCode())
        {
            case 10:
            case 20:
            case 30:
            case 40:
            case 210:
            case 220:
            case 230:
            {
                std::istringstream iss(isit->GetValue());
                iss >> value;
                if (iss.fail())
                {
                    throw DxfFileParserException("Error while reading value.");
                }
                break;
            }

            default:
                break;
        }

        switch (isit->GetCode())
        {
            case 10:  centreX = value; break;
            case 20:  centreY = value; break;
            case 30:  centreZ = value; break;
            case 40:  radius  = value; break;
            case 210: normalX = value; break;
            case 220: normalY = value; break;
            case 230: normalZ = value; break;
            default: break;
        }

        ++isit;
    }

    return std::make_shared<DxfEntityCircle>
    (
        point(centreX, centreY, centreZ),
        radius,
        vector(normalX, normalY, normalZ)
    );
}


std::shared_ptr<DxfEntity> DxfFileParser::ParseEllipse
(
    DxfFileParser::IsItType& isit
)
{
    double centreX = 0;
    double centreY = 0;
    double centreZ = 0;
    double endMajorX = 0;
    double endMajorY = 0;
    double endMajorZ = 0;
    double normalX = 0;
    double normalY = 0;
    double normalZ = 1;
    double axisRatio  = 0;
    double startAngle = 0;
    double endAngle   = 0;

    ++isit;
    while (isit->GetCode() != 0)
    {
        double value;

        switch (isit->GetCode())
        {
            case 10:
            case 20:
            case 30:
            case 11:
            case 21:
            case 31:
            case 210:
            case 220:
            case 230:
            case 40:
            case 41:
            case 42:
            {
                std::istringstream iss(isit->GetValue());
                iss >> value;
                if (iss.fail())
                {
                    throw DxfFileParserException("Error while reading value.");
                }
                break;
            }

            default:
                break;
        }

        switch (isit->GetCode())
        {
            case 10:  centreX = value;    break;
            case 20:  centreY = value;    break;
            case 30:  centreZ = value;    break;
            case 11:  endMajorX = value;  break;
            case 21:  endMajorY = value;  break;
            case 31:  endMajorZ = value;  break;
            case 210: normalX = value;    break;
            case 220: normalY = value;    break;
            case 230: normalZ = value;    break;
            case 40:  axisRatio = value;  break;
            case 41:  startAngle = value; break;
            case 42:  endAngle = value;   break;
            default: break;
        }

        ++isit;
    }

    return std::make_shared<DxfEntityEllipse>
    (
        point(centreX, centreY, centreZ),
        point(endMajorX, endMajorY, endMajorZ),
        vector(normalX, normalY, normalZ),
        axisRatio,
        startAngle,
        endAngle
    );
}


std::shared_ptr<DxfEntity> DxfFileParser::ParseInsert
(
    DxfFileParser::IsItType& isit
)
{
    std::string blockName;

    double xInsertionPoint = 0;
    double yInsertionPoint = 0;
    double zInsertionPoint = 0;

    double xScaleFactor = 1;
    double yScaleFactor = 1;
    double zScaleFactor = 1;

    double rotationAngle = 0;

    int columnCount = 1;
    int rowCount    = 1;

    double columnSpacing = 0;
    double rowSpacing    = 0;

    double xExtrusionDirection = 0;
    double yExtrusionDirection = 0;
    double zExtrusionDirection = 1;

    ++isit;
    while (isit->GetCode() != 0)
    {
        std::string valueAsString;
        double      valueAsDouble;
        int         valueAsInt;

        switch (isit->GetCode())
        {
            case 2:
            {
                std::istringstream iss(isit->GetValue());
                iss >> valueAsString;
                if (iss.fail())
                {
                    throw DxfFileParserException("Error while reading value.");
                }
                break;
            }

            case 10:
            case 20:
            case 30:
            case 41:
            case 42:
            case 43:
            case 50:
            case 44:
            case 45:
            case 210:
            case 220:
            case 230:
            {
                std::istringstream iss(isit->GetValue());
                iss >> valueAsDouble;
                if (iss.fail())
                {
                    throw DxfFileParserException("Error while reading value.");
                }
                break;
            }

            case 70:
            case 71:
            {
                std::istringstream iss(isit->GetValue());
                iss >> valueAsInt;
                if (iss.fail())
                {
                    throw DxfFileParserException("Error while reading value.");
                }
                break;
            }
        }

        switch (isit->GetCode())
        {
            case 2:   blockName           = valueAsString; break;
            case 10:  xInsertionPoint     = valueAsDouble; break;
            case 20:  yInsertionPoint     = valueAsDouble; break;
            case 30:  zInsertionPoint     = valueAsDouble; break;
            case 41:  xScaleFactor        = valueAsDouble; break;
            case 42:  yScaleFactor        = valueAsDouble; break;
            case 43:  zScaleFactor        = valueAsDouble; break;
            case 50:  rotationAngle       = valueAsDouble; break;
            case 44:  columnSpacing       = valueAsDouble; break;
            case 45:  rowSpacing          = valueAsDouble; break;
            case 210: xExtrusionDirection = valueAsDouble; break;
            case 220: yExtrusionDirection = valueAsDouble; break;
            case 230: zExtrusionDirection = valueAsDouble; break;
            case 70:  columnCount         = valueAsInt;    break;
            case 71:  rowCount            = valueAsInt;    break;
        }

        ++isit;
    }

    return std::make_shared<DxfEntityInsert>
    (
        blockName,
        xInsertionPoint,
        yInsertionPoint,
        zInsertionPoint,
        xScaleFactor,
        yScaleFactor,
        zScaleFactor,
        rotationAngle,
        columnCount,
        rowCount,
        columnSpacing,
        rowSpacing,
        xExtrusionDirection,
        yExtrusionDirection,
        zExtrusionDirection
    );
}


std::shared_ptr<DxfEntity> DxfFileParser::ParseLine
(
    DxfFileParser::IsItType& isit
)
{
    double startX  = 0;
    double startY  = 0;
    double startZ  = 0;
    double endX    = 0;
    double endY    = 0;
    double endZ    = 0;
    double normalX = 0;
    double normalY = 0;
    double normalZ = 1;

    ++isit;
    while (isit->GetCode() != 0)
    {
        double value;

        switch (isit->GetCode())
        {
            case 10:
            case 20:
            case 30:
            case 11:
            case 21:
            case 31:
            case 210:
            case 220:
            case 230:
            {
                std::istringstream iss(isit->GetValue());
                iss >> value;
                if (iss.fail())
                {
                    throw DxfFileParserException("Error while reading value.");
                }
                break;
            }

            default:
                break;
        }

        switch (isit->GetCode())
        {
            case 10:  startX = value;  break;
            case 20:  startY = value;  break;
            case 30:  startZ = value;  break;
            case 11:  endX = value;    break;
            case 21:  endY = value;    break;
            case 31:  endZ = value;    break;
            case 210: normalX = value; break;
            case 220: normalY = value; break;
            case 230: normalZ = value; break;
            default: break;
        }

        ++isit;
    }

    return std::make_shared<DxfEntityLine>
    (
        point(startX, startY, startZ),
        point(endX, endY, endZ),
        vector(normalX, normalY, normalZ)
    );
}


std::shared_ptr<DxfEntity> DxfFileParser::ParseLwPolyLine
(
    DxfFileParser::IsItType& isit
)
{
    int numOfVertices = 0;
    int polyLineFlag = 0;
    std::vector<double> verticesX;
    std::vector<double> verticesY;
    double normalX = 0;
    double normalY = 0;
    double normalZ = 1;

    double elevation = 0.0;

    ++isit;
    while (isit->GetCode() != 0)
    {
        double value;

        switch (isit->GetCode())
        {
            case 90:
            case 70:
            case 10:
            case 20:
            case 210:
            case 220:
            case 230:
            {
                std::istringstream iss(isit->GetValue());
                iss >> value;
                if (iss.fail())
                {
                    throw DxfFileParserException("Error while reading value.");
                }
                break;
            }

            default:
                break;
        }

        switch (isit->GetCode())
        {
            case 90:  numOfVertices = value;      break;
            case 70:  polyLineFlag = value;       break;
            case 10:  verticesX.push_back(value); break;
            case 20:  verticesY.push_back(value); break;
            case 210: normalX = value;            break;
            case 220: normalY = value;            break;
            case 230: normalZ = value;            break;
            case 38:  elevation = value;          break;
            default: break;
        }

        ++isit;
    }

    if
    (
        static_cast<int>(verticesX.size()) != numOfVertices
     || static_cast<int>(verticesY.size()) != numOfVertices)
    {
        throw DxfFileParserException
        (
            "Invalid number of parsed points of LWPOLYLINE entity."
        );
    }

    return std::make_shared<DxfEntityPolyLine>
    (
        numOfVertices,
        polyLineFlag,
        verticesX,
        verticesY,
        vector(normalX, normalY, normalZ),
        elevation
    );
}


std::shared_ptr<DxfEntity> DxfFileParser::ParseSpline
(
    DxfFileParser::IsItType& isit
)
{
    double normalX = 0;
    double normalY = 0;
    double normalZ = 0;

    int    splineFlag = 0;
    int    degree = 0;
    int    numOfKnots = 0;
    int    numOfControlPoints = 0;
    int    numOfFitPoints = 0;
    double knotTol = 0.0000001;
    double controlPointTol = 0.0000001;
    double fitTol = 0.0000000001;

    double startTangentX = 0;
    double startTangentY = 0;
    double startTangentZ = 0;

    double endTangentX = 0;
    double endTangentY = 0;
    double endTangentZ = 0;

    std::list<double> knotValues;
    std::list<double> weights;

    std::list<double> controlPointsX;
    std::list<double> controlPointsY;
    std::list<double> controlPointsZ;

    std::list<double> fitPointsX;
    std::list<double> fitPointsY;
    std::list<double> fitPointsZ;

    ++isit;
    while (isit->GetCode() != 0)
    {
        double value;

        switch (isit->GetCode())
        {
            case 210:
            case 220:
            case 230:
            case 70:
            case 71:
            case 72:
            case 73:
            case 74:
            case 42:
            case 43:
            case 44:
            case 12:
            case 22:
            case 32:
            case 13:
            case 23:
            case 33:
            case 40:
            case 41:
            case 10:
            case 20:
            case 30:
            case 11:
            case 21:
            case 31:
            {
                std::istringstream iss(isit->GetValue());
                iss >> value;
                if (iss.fail())
                {
                    throw DxfFileParserException("Error while reading value.");
                }
                break;
            }

            default:
                break;
        }

        switch (isit->GetCode())
        {
            case 210: normalX = value;                  break;
            case 220: normalY = value;                  break;
            case 230: normalZ = value;                  break;
            case 70:  splineFlag = value;               break;
            case 71:  degree = value;                   break;
            case 72:  numOfKnots = value;               break;
            case 73:  numOfControlPoints = value;       break;
            case 74:  numOfFitPoints = value;           break;
            case 42:  knotTol = value;                  break;
            case 43:  controlPointTol = value;          break;
            case 44:  fitTol = value;                   break;
            case 12:  startTangentX = value;            break;
            case 22:  startTangentY = value;            break;
            case 32:  startTangentZ = value;            break;
            case 13:  endTangentX = value;              break;
            case 23:  endTangentY = value;              break;
            case 33:  endTangentZ = value;              break;
            case 40:  knotValues.push_back(value);      break;
            case 41:  weights.push_back(value);         break;
            case 10:  controlPointsX.push_back(value);  break;
            case 20:  controlPointsY.push_back(value);  break;
            case 30:  controlPointsZ.push_back(value);  break;
            case 11:  fitPointsX.push_back(value);      break;
            case 21:  fitPointsY.push_back(value);      break;
            case 31:  fitPointsZ.push_back(value);      break;
            default: break;
        }

        ++isit;
    }

    if (static_cast<int>(knotValues.size()) != numOfKnots)
    {
        throw DxfFileParserException
        (
            "Invalid number of knots of LWPOLYLINE entity."
        );
    }

    // if (weights.size() != numOfKnots) <<< ?????
    // {
    //     throw DxfFileParserException
    //     (
    //         "Invalid number of weights of LWPOLYLINE entity."
    //     );
    // }

    if
    (
        static_cast<int>(controlPointsX.size()) != numOfControlPoints
     || static_cast<int>(controlPointsY.size()) != numOfControlPoints
     || static_cast<int>(controlPointsZ.size()) != numOfControlPoints
    )
    {
        throw DxfFileParserException
        (
            "Invalid number of control points of LWPOLYLINE entity."
        );
    }

    if
    (
        static_cast<int>(fitPointsX.size()) != numOfFitPoints
     || static_cast<int>(fitPointsY.size()) != numOfFitPoints
     || static_cast<int>(fitPointsZ.size()) != numOfFitPoints
    )
    {
        throw DxfFileParserException
        (
            "Invalid number of fit points of LWPOLYLINE entity."
        );
    }

    return std::make_shared<DxfEntitySpline>
    (
        vector(normalX, normalY, normalZ),
        splineFlag,
        degree,
        numOfKnots,
        numOfControlPoints,
        numOfFitPoints,
        knotTol,
        controlPointTol,
        fitTol,
        vector(startTangentX, startTangentY, startTangentZ),
        vector(endTangentX, endTangentY, endTangentZ),
        knotValues,
        weights,
        controlPointsX,
        controlPointsY,
        controlPointsZ,
        fitPointsX,
        fitPointsY,
        fitPointsZ
    );
}

} // End namespace Foam
