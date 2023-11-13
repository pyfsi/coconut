#ifndef DxfFileParser_H
#define DxfFileParser_H

#include "DxfFile.H"

#include <istream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>

namespace Foam
{

// ============================================================================
// class DxfFileParserException
// ============================================================================

class DxfFileParserException : public std::runtime_error
{
    using std::runtime_error::runtime_error;
};


// ============================================================================
// class DxfRecord
// ============================================================================

class DxfRecord
{
public:

    // ------------------------------------------------------------------------
    // Public member functions
    // ------------------------------------------------------------------------

    int GetCode() const
    {
        return code_;
    }

    void SetCode(int code)
    {
        code_ = code;
    }

    const std::string& GetValue() const
    {
        return value_;
    }

    void SetValue(const std::string& value)
    {
        value_ = value;
    }

private:

    // ------------------------------------------------------------------------
    // Private data members
    // ------------------------------------------------------------------------

    int         code_;
    std::string value_;

};


// ============================================================================
// Non-member functions
// ============================================================================

std::istream& operator>>(std::istream& is, DxfRecord& record);


// ============================================================================
// class DxfFileParser
// ============================================================================

class DxfFileParser
{
public:

    // ------------------------------------------------------------------------
    // Public member functions
    // ------------------------------------------------------------------------

    void Parse(std::istream& is, DxfFile& file);

private:

    // ------------------------------------------------------------------------
    // Private member types
    // ------------------------------------------------------------------------

    using IsItType = std::istream_iterator<DxfRecord>;

    // ------------------------------------------------------------------------
    // Private member functions
    // ------------------------------------------------------------------------

    //
    // Section HEADER
    //

    void ParseSectionHeader(IsItType& isit, DxfFile& file);

    DxfFile::Units      ParseUnits(IsItType& isit);
    DxfFile::AngleUnits ParseAngleUnits(IsItType& isit);

    //
    // Section BLOCKS
    //

    void ParseSectionBlocks(IsItType& isit, DxfFile& file);

    std::shared_ptr<DxfBlock> ParseBlock(IsItType& isit);

    //
    // Section ENTITIES
    //

    void ParseSectionEntities(IsItType& isit, DxfFile& file);

    std::shared_ptr<DxfEntity> ParseArc(IsItType& isit);
    std::shared_ptr<DxfEntity> ParseCircle(IsItType& isit);
    std::shared_ptr<DxfEntity> ParseEllipse(IsItType& isit);
    std::shared_ptr<DxfEntity> ParseInsert(IsItType& isit);
    std::shared_ptr<DxfEntity> ParseLine(IsItType& isit);
    std::shared_ptr<DxfEntity> ParseLwPolyLine(IsItType& isit);
    std::shared_ptr<DxfEntity> ParseSpline(IsItType& isit);

};

} // End namespace Foam

#endif // DxfFileParser_H
