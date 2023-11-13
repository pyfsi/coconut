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
    rollingMillMesh

Description
    Creates 3D meshes for rolling mill components

SourceFiles
    rollingMillMesh.C

\*---------------------------------------------------------------------------*/

#ifndef rollingMillMesh_H
#define rollingMillMesh_H

#include "polyMeshGen.H"
#include "IOdictionary.H"
#include "workflowControls.H"
#include "triSurf.H"

#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class meshSurfaceEngine;
class rollerSurfaceCreator;
class Time;
class rollingMillContactHandler;
class rollingMillGeometryHandler;
class rollingMillPatchNamesHandler;

/*---------------------------------------------------------------------------*\
                    Class dieGeometryInfo Declaration
\*---------------------------------------------------------------------------*/

class dieGeometryInfo
{
    // Private data
        //- pointer to the dictionary for the die
        const dictionary* dieDictPtr_;

        //- axial length of a die
        scalar axialLength_;

        //- outer radius of a die
        scalar outerRadius_;

        //- inner radius at the inlet
        scalar inletRadius_;

        //- outer radius at the inlet
        //- it is not used if the surface mesh of the profile is available
        scalar outletRadius_;

        //- geometry tolerance
        scalar geometryTolerance_;

        //- wedge angle used for axi-symmetry
        scalar wedgeAngle_;

        //- surface mesh of the axial profile of a die
        const triSurf* innerSurfPtr_;

        //- cross-sections at various axial positions
        PtrList<std::pair<scalar, std::shared_ptr<triSurf> > > crossSections_;

        //- surface mesh of the inlet profile
        //- positioned in the y-z plane
        const triSurf* inletSurfPtr_;

        //- surface mesh of the outlet profile
        //- positioned in the y-z plane
        const triSurf* outletSurfPtr_;

        //- surface mesh of the casing
        const triSurf* casingSurfPtr_;

        //- const reference to the rolling mill geometry handler
        const rollingMillGeometryHandler& geomHandler_;

        //- reference to the patch handler
        const rollingMillPatchNamesHandler& patchHandler_;

public:

    // Constructor

        //- Construct from components
        dieGeometryInfo
        (
            const dictionary& dieDict,
            const rollingMillGeometryHandler& geomHandler,
            const rollingMillPatchNamesHandler& patchHandler_
        );

    // Destructor

        ~dieGeometryInfo();

    // Public member functions
        //- reference to the die dictionary
        const dictionary& dieDict() const;

        //- return the outer diameter of a die
        scalar outerDiameter() const;

        //- return the average diameter and the wire inlet
        scalar inletDiameter() const;

        //- is the inlet circular or not
        bool isInletCircular() const;

        //- is the outlet circular or not
        bool isOutletCircular() const;

        //- return theaverage diameter at the outlet
        scalar outletDiameter() const;

        //- is it a profiled die or not
        bool isProfiled() const;

        //- is the setup symmetric or not
        bool isSymmetric() const;

        //- is casing present or not
        bool isCasingPresent() const;

        //- return type of symmetry
        direction typeOfSymmetry() const;

        //- is it a wedge or not
        bool isWedge() const;

        //- return wedge angle
        scalar wedgeAngle() const;

        //- start circumferential angle
        scalar startCircumAngle() const;

        //- end circumferential angle
        scalar endCircumAngle() const;

        //- const reference to the surface mesh for profiled dies at inlet
        const triSurf& inletSurface() const;

        //- const reference to the surface mesh for profiled dies at outlet
        const triSurf& outletSurface() const;

        //- const reference to user-provided cross section positions
        //- each pair represents a surface mesh at a given axial coordinate
        //- given in absolute values
        const PtrList<std::pair<scalar, std::shared_ptr<triSurf> > >&
        crossSections() const;

        //- axial cross-section of the geometry in the x-y plane
        //- used for generating a 2D mesh of a die
        const triSurf *axialCrossSectionSurface() const;

        //- axial cross-section of the geometry in the x-y plane
        //- used for generating a 2D mesh of a casing
        const triSurf* casingCrossSectionSurface() const;
};

/*---------------------------------------------------------------------------*\
                    Class wireGeometryInfo Declaration
\*---------------------------------------------------------------------------*/

class wireGeometryInfo
{
    // Private data
        //- surface mesh positioned correctly in the y-z plane
        const triSurf* surfPtr_;

        //- pointer to wireMeshDict
        const dictionary* wireDictPtr_;

        //- outer diameter of the wire
        scalar wireDiameter_;

        //- axial grading
        scalar axialGrading_;

        //- axial shift
        vector axialShift_;

        //- resolution in the axial direction
        label axialResolution_;

        //- is it a rolling pass
        bool isRollingPass_;

        //- is it a drawing pass
        bool isDrawingPass_;

        //- name of the wall patch
        word wallPatchName_;

        //- type of symmetry determined by wither die or a roller setup
        direction symmetryType_;

        //- wedge angle for axi-symmetry
        scalar wedgeAngle_;

        //- shall the coating be meshed or not
        bool isCoatingPresent_;

public:

    // Static data
        //- transformation matrix into the x-y coordinates for meshing
        static const tensor transformationMatrix;

    // Constructor

        //- Construct from components
        wireGeometryInfo
        (
            const dictionary& wireDict,
            const triSurf* surfacePtr = NULL
        );

    // Destructor

        ~wireGeometryInfo();

    // Public member functions
        //- is cross-section circular
        bool isCrossSectionCircular() const;

        //- const reference to the surface mesh
        const triSurf& surface() const;

        //- reference to the wire dictionary
        const dictionary& wireDict() const;

        //- get and set wire diameter
        scalar& wireDiameter();
        const scalar& wireDiameter() const;

        //- get and set axial grading
        const scalar& axialGrading() const;
        scalar& axialGrading();

        //- get and set axial resolution
        label axialResolution() const;
        label& axialResolution();

        //- axial shift
        const vector& axialShift() const;
        vector& axialShift();

        //- get and set the name of th wall patch
        const word& wallPatchName() const;
        word& wallPatchName();

        //- get and set the type of symmetry
        direction typeOfSymmetry() const;
        direction& typeOfSymmetry();

        //- wedge angle for axi-symmetry
        scalar wedgeAngle() const;
        scalar& wedgeAngle();

        //- is coating present on the wire or not
        bool isCoatingPresent() const;
        bool& isCoatingPresent();

        //- a pointer to the transformed surface mesh used for
        //- 2D meshing in the x-y plane
        const triSurf* transformedSurface() const;

        //- transform the mesh into requested position
        void backTransform2DMeshToPosition(polyMeshGen& mesh) const;
};

/*---------------------------------------------------------------------------*\
                    Class rollGeometryInfo Declaration
\*---------------------------------------------------------------------------*/

class rollGeometryInfo
{
    // Private data
        //- surface mesh positioned correctly in the y-z plane
        mutable triSurf* surfPtr_;

        //- pointer to the roll dictionary
        const dictionary* rollDictPtr_;

        //- outer diameter
        scalar outerDiameter_;

        //- roll position
        word rollPosition_;

        //- roll offset
        word rollOffset_;

        //- grading in the circumferential direction
        scalar circumGrading_;

        //- scaling in the circumferential direction
        scalar circumScaling_;

        //- number of subdivisions in the circumferential direction
        label circumferentialResolution_;

        //- number of subdivisions in the axial direction
        label axialResolution_;

        //- translation in the axial direction
        scalar axialTranslation_;

        //- translation in the radial direction
        scalar radialTranslation_;

        //- grading factor in the y direction
        scalar radialGrading_;

        //- grading factor in the z direction
        scalar axialGrading_;

        //- scaling factor in the axial direction
        scalar axialScaling_;

        //- scaling factor in the radial direction
        scalar radialScaling_;

        //- translation in the x-direction
        scalar offsetX_;

        //- translation in the radial direction
        scalar radialOffset_;

        //- width of the contact area centered at origin
        scalar contactWidth_;

        //- extra angle added to the contact patch
        scalar extraContactAngle_;

        //- symmetry type of a roll
        direction symmetryType_;

        //- closed roll setup
        bool isClosedSetup_;

        //- centre roll in axial direction
        bool centreInAxialDirection_;

    // Static private member functions
        //- populate rollPositions with values
        static wordHashSet populateRollPositions();

public:

    // Static data
        //- available roll positions
        static const wordHashSet rollPositions;

        //- transformation matrix into the x-y coordinates for meshing
        static const tensor transformationMatrix;

    // Constructor

        //- Construct from components
        rollGeometryInfo
        (
            triSurf* surfacePtr,
            const dictionary& rollDict,
            const word rollPosition
        );

    // Destructor

        ~rollGeometryInfo();

    // Public member functions
        //- const reference to the surface mesh
        const triSurf& surface() const;

        //- const reference to the dictionary for the given roll
        const dictionary& rollDict() const;

        //- roll offset name
        word& rollOffset();
        const word& rollOffset() const;

        //- get and set outer diameter
        scalar& outerDiameter();
        const scalar& outerDiameter() const;

        //- get and set the grading in the circumferential direction
        scalar& circumGrading();
        const scalar& circumGrading() const;

        //- scaling factor in the circumferential direction
        scalar& circumScaling();
        const scalar& circumScaling() const;

        //- get and set circumferential resolution
        label& circumResolution();
        label circumResolution() const;

        //- get and set axial resolution
        label& axialResolution();
        label axialResolution() const;

        //- get and set axial shift distance
        scalar& axialShiftDistance();
        const scalar& axialShiftDistance() const;

        //- get and set radial shift distance
        scalar& radialShiftDistance();
        const scalar& radialShiftDistance() const;

        //- grading in the y diretion
        scalar& radialGrading();
        const scalar& radialGrading() const;

        //- grading in the z direction
        scalar& axialGrading();
        const scalar& axialGrading() const;

        //- radial scaling
        scalar& radialScaling();
        const scalar& radialScaling() const;

        //- axial scaling
        scalar& axialScaling();
        const scalar& axialScaling() const;

        //- translation of the roll in the x-direction
        scalar& offsetX();
        const scalar& offsetX() const;

        //- translation of the roll in the radial direction
        scalar& radialOffset();
        const scalar& radialOffset() const;

        //- get the range of the contact patch specified by the user
        //- it is negative if it is not set
        scalar& contactWidth();
        const scalar& contactWidth() const;

        //- extra angle for the contact patch
        scalar& extraContactAngle();
        scalar extraContactAngle() const;

        //- get and set symmetry configuration
        direction& typeOfSymmetry();
        direction typeOfSymmetry() const;

        //- is it closed setup (rolls touch each other)
        bool& isClosedSetup();
        bool isClosedSetup() const;

        //- shall the roll be centered in the axial direction
        bool& centreInAxialDirection();
        bool centreInAxialDirection() const;

        //- a pointer to the transformed surface mesh used for
        //- 2D meshing in the x-y plane
        const triSurf* transformedSurface() const;

        //- transform the mesh into requested position
        void backTransform2DMeshToPosition(polyMeshGen& mesh) const;

        //- roll position
        word rollPosition() const;

        //- rotation axis
        vector rotationAxis() const;

        //- translation distance in the axial direction
        vector axialShift() const;

        //- a vector in the radial direction
        //- the whole roller shall be shifted in that direction
        //- depends on the position of the roller
        vector radialShift(const scalar wireDiameter = 0.0) const;

        //- a point on the rotation axis
        point origin() const;

        //- grading and scaling requirements
        dictionary anisotropic2DMeshing() const;

    // Modifiers
        //- update surface patches
        void updateSurfacePatches
        (
            const wordList& patchNames,
            const wordList& patchTypes,
            const labelLongList& patchIndices
        ) const;
};

/*---------------------------------------------------------------------------*\
                Class rollingMillGeometryHandler Declaration
\*---------------------------------------------------------------------------*/

class rollingMillGeometryHandler
{
    // Private data
        //- reference to runTime
        const Time& runTime_;

        //- reference to the patch handler
        const rollingMillPatchNamesHandler& patchHandler_;

        //- a list of roll positions
        PtrList<rollGeometryInfo> geometries_;

        //- wire geometry if needed
        wireGeometryInfo* wireGeomPtr_;

        //- die geometry if needed
        dieGeometryInfo* dieGeomPtr_;

        //- wire length
        scalar wireLength_;

        //- wire diameter
        scalar wireDiameter_;

        //- contact length
        scalar contactLength_;

        //- extra angle for the contact patch
        scalar extraContactAngle_;

        //- initial clearance between rolls and origin
        scalar initialRollGapClearance_;

        //- maximum allowed deviation of the discrete geometry from
        //- the analytical one
        scalar geometryTolerance_;

        //- is it the roller setup
        bool isRollerSetup_;

        //- it is a closed roll setup
        bool areRollsClosed_;

        //- position the centre of rolls in their axial direction
        bool centreInAxialDirection_;

        //- generate only one contact patch
        bool singleRollerContactPatch_;

        //- type of symmetry
        direction symmetryType_;

    // Private member functions
        //- check if all points are located in a positive or negative semi space
        bool checkHalfSymmetry() const;

        //- check if all points are locates in a single quadrant of the yz plane
        bool checkQuarterSymmetry() const;

        //- write a dictionary with roll positions in operation mode
        void writeRollPositions() const;

        //- position the roll given by the point coordinates
        //- into the correct position
        void transformRollToPosition
        (
            triSurf& rollSurf,
            const word& rollPosition,
            const bool flipRoll
        ) const;

        //- cut the geometry to enforce symmetry as specified by the user
        bool cutGeometry(triSurf& rollSurf, const word& rollPosition) const;

        //- parse roll dictionary and generate the surface for meshing
        void parseRollDictionary
        (
            const dictionary& rollDict,
            const word rollOffsetName = ""
        );

        //- parse dictionaries
        void parseDictionaries(const dictionary& meshDict);

        //- disallow bitwise asignment
        void operator=(const rollingMillGeometryHandler&);

        //- disallow copy construct
        rollingMillGeometryHandler(const rollingMillGeometryHandler&);

public:

    // Static data
        //- allowed types of symmetric configurations
        static const wordHashSet symmetryTypes;

    // Static member functions

        //- populate names of symmetric configurations
        static wordHashSet populateSymmetricConfigurations();

    // Enumerators
        //- types of symmetry
        enum symmetryTypes_
        {
            NONE = 0,
            POSY = 1,
            NEGY = 2,
            POSZ = 4,
            NEGZ = 8,
            AXISYMMETRIC = 16
        };

    // Constructors
        //- construct from dictionary
        rollingMillGeometryHandler
        (
            const Time& runTime,
            const dictionary& meshDict,
            const rollingMillPatchNamesHandler& patchHandler
        );

    // Destructor

        ~rollingMillGeometryHandler();

    // Public member functions
        //- check the existence of symmetry in the roll setup
        direction symmetryType() const;

        //- return roll positions
        const PtrList<rollGeometryInfo>& rollPositions() const;

        //- return the wire geometry info
        const wireGeometryInfo& wireGeometry() const;

        //- return the die geometry info
        const dieGeometryInfo& dieGeometry() const;

        //- return geometry tolerance
        scalar geometryTolerance() const;

        //- return wire length
        scalar wireLength() const;

        //- return wire diameter
        scalar wireDiameter() const;

        //- return contact length
        scalar contactLength() const;

        //- extra angle for the contact patch
        scalar extraAngle() const;

        //- initial clearance between rolls and origin
        scalar initialRollGapClearance() const;

        //- are rolls set up or not
        bool isRollSetupValid() const;

        //- is it a closed setup of rolls or not (rolls touch each other)
        bool areRollsClosed() const;

        //- set to true if only one contact patch shall be generated
        bool singleRollerContactPatch() const;

        //- parse wire dictionary
        void parseWireDictionary(const dictionary& meshDict);

        //- parse die dictionary
        void parseDieDictionary(const dictionary& meshDict);
};

/*---------------------------------------------------------------------------*\
                Class rollingMillPatchNamesHandler Declaration
\*---------------------------------------------------------------------------*/

class rollingMillPatchNamesHandler
{
    // Private data
        //- const reference to meshDict
        const dictionary& meshDict_;

        //- roller patch names and types
        std::map<label, word> rollerPatchNames_;

        //- wire patch names
        std::map<label, word> wirePatchNames_;

        //- wire coating patch names
        std::map<label, word> coatingPatchNames_;

        //- die patch names
        std::map<label, word> diePatchNames_;

        //- casing patch names
        std::map<label, word> casingPatchNames_;

    // Private member functions
        //- default roller patch names
        void defaultRollerPatchNames();

        //- default patch names for wire meshes
        void defaultWirePatchNames();

        //- default patch names for wire coating
        void defaultWireCoatingPatchNames();

        //- default patch names for die meshes
        void defaultDiePatchNames();

        //- default patch names for casing meshes
        void defaultCasingPatchNames();

        //- disallow bitwise asignment
        void operator=(const rollingMillPatchNamesHandler&);

        //- disallow copy construct
        rollingMillPatchNamesHandler(const rollingMillPatchNamesHandler&);

public:

    // Enumerators
        //- patch keys for rollers
        enum rollerPatchKeys
        {
            ROLLERCONTACT = 0,
            ROLLERAXIS,
            ROLLERFRONT,
            ROLLERBACK,
            ROLLERSYMMY,
            ROLLERSYMMZ,
            ROLLERTOROLLERBACK,
            ROLLERTOROLLERFRONT,
            ROLLERTOAIR
        };

        //- patch keys for wire meshes
        enum wirePatchKeys
        {
            WIRECONTACT = 0,
            WIRECOATING,
            WIREUPSTREAM,
            WIREDOWNSTREAM,
            WIRESYMMY,
            WIRESYMMZ,
            WIREFRONT,
            WIREBACK
        };

        //- patch keys for wire coating
        enum coatingPatchKeys
        {
            COATINGCONTACT = 0,
            COATINGWIRE,
            COATINGUPSTREAM,
            COATINGDOWNSTREAM,
            COATINGSYMMY,
            COATINGSYMMZ,
            COATINGFRONT,
            COATINGBACK
        };

        //- patch keys for die meshes
        enum diePatchKeys
        {
            DIEUPSTREAM = 0,
            DIEDOWNSTREAM,
            DIEWIRE,
            DIEHOUSING,
            DIESYMMY,
            DIESYMMZ,
            DIEFRONT,
            DIEBACK,
            DIEENTRANCECONE,
            DIEEXITCONE
        };

        //- patch keys for casing meshes
        enum casingPatchKeys
        {
            CASINGDOWNSTREAM = 0,
            CASINGUPSTREAM,
            CASINGTOOUTSIDE,
            CASINGENTRANCE,
            CASINGTODIERADIAL,
            CASINGTODIEAXIAL,
            CASINGEXIT,
            CASINGSYMMY,
            CASINGSYMMZ,
            CASINGFRONT,
            CASINGBACK
        };

    // Constructors
        //- construct from dictionary
        rollingMillPatchNamesHandler(const dictionary& meshDict);

    // Destructor

        ~rollingMillPatchNamesHandler();

    // Public member functions
        //- return patch name for roll
        word patchNameForRoll(const rollerPatchKeys&) const;

        //- return patch type for roll
        word patchTypeForRoll
        (
            const rollerPatchKeys&,
            const rollingMillGeometryHandler::symmetryTypes_&
        ) const;

        //- return patch name for wire
        word patchNameForWire(const wirePatchKeys&) const;

        //- return patch type for wire
        word patchTypeForWire
        (
            const wirePatchKeys&,
            const rollingMillGeometryHandler::symmetryTypes_&
        ) const;

        //- return patch name for wire coating
        word patchNameForCoating(const coatingPatchKeys&) const;

        //- return patch type for wire coating
        word patchTypeForCoating
        (
            const coatingPatchKeys&,
            const rollingMillGeometryHandler::symmetryTypes_&
        ) const;

        //- return patch name for die
        word patchNameForDie(const diePatchKeys&) const;

        //- return patch type for die
        word patchTypeForDie
        (
            const diePatchKeys&,
            const rollingMillGeometryHandler::symmetryTypes_&
        ) const;

        //- return patch name for casing
        word patchNameForCasing(const casingPatchKeys&) const;

        //- return patch type for casing
        word patchTypeForCasing
        (
            const casingPatchKeys&,
            const rollingMillGeometryHandler::symmetryTypes_&
        ) const;
};

/*---------------------------------------------------------------------------*\
                Class rollingMillContactHandler Declaration
\*---------------------------------------------------------------------------*/

class rollingMillContactHandler
{
    // Private data
        //- const reference to a dictionary
        const dictionary& meshDict_;

        //- const reference to roll geometries
        const PtrList<rollGeometryInfo>& rollGeometries_;

        //- const reference to patch names handler
        const rollingMillPatchNamesHandler& patchNamesHandler_;

        //- write cross section of the rolls
        Switch writeRollProfile_;

        //- distance tolerance for evaluating contacts
        scalar contactDistanceTolerance_;

        //- cross section edges
        edgeLongList crossSectionEdges_;

        //- point label mapping from the volume mesh onto the cross section
        labelLongList pointMapping_;

        //- surface mesh of the cross-section
        triSurf* crossSectionSurfPtr_;

        //- manifold for the surface triangle
        //- each roll is a separate manifold
        labelLongList manifold_;

        //- roll position to manifold mapping
        std::map<word, label> rollPositionToManifold_;

        //- octree for searches
        meshOctree* octreePtr_;

    // Private member functions
        //- parse dictionary and get the roll gaps
        bool readGeometryData();

        //- create the octree of the cross-section
        void createOctree();

        //- position cross sections
        void positionCrossSections();

        //- position dependent cross sections
        //- this hapens in case of side rolls setup
        void positionDependentCrossSections();

        //- search for geometry parts within close contact to each other
        void checkForContactRegions();

        //- disallow bitwise asignment
        void operator=(const rollingMillContactHandler&);

        //- disallow copy construct
        rollingMillContactHandler(const rollingMillContactHandler&);

    // Enumerators

        enum patchType
        {
            CONTACT = 1,
            AXIS = 2,
            SYMMETRY = 4,
            FRONT = 8,
            BACK = 16
        };

        enum regionType
        {
            NONE = 0,
            CONTACTBACK = 1,
            CONTACTFRONT = 2,
            ROLLERTOAIR = 4,
            ROLLERTOROLLERCONTACTPATCH = 8,
            ROLLERTOROLLERCONTACTFRONT = 16,
            ROLLERTOROLLERCONTACTBACK = 32,
            ROLLERTOWIRECONTACTPATCH = 64
        };

public:

    // Constructors
        //- construct from mesh and dictionary
        rollingMillContactHandler
        (
            const dictionary& meshDict,
            const PtrList<rollGeometryInfo>& rollGeometries,
            const rollingMillPatchNamesHandler&
        );

    // Destructor

        ~rollingMillContactHandler();

    // Public member functions
};

/*---------------------------------------------------------------------------*\
                    Class rollingMillMesh Declaration
\*---------------------------------------------------------------------------*/

class rollingMillMesh
{
    // Private data
        //- reference to Time
        const Time& db_;

        //- region name
        const fileName regionName_;

        //- time step
        const fileName timeStep_;

        //- IOdictionary containing information about cell sizes, etc..
        IOdictionary meshDict_;

        //- global mesh of the whole setup
        polyMeshGen* globalMeshPtr_;

        //- patch names handler
        rollingMillPatchNamesHandler patchHandler_;

        //- geometry handler
        rollingMillGeometryHandler geomHandler_;

    // Private member functions
        //- generate a 3D wire mesh based on the data in the dictionary
        void generateWireMesh();

        //- generate a 3D mesh of the wire coating
        void generateWireCoatingMesh();

        //- generate 3D mesh from a cross-section mesh
        void generateRollerMeshes();

        //- generate 3D mesh of a die
        void generateDieMesh();

        //- generate 3D mesh of a casing
        void generateCasingMesh();

    // Private nested classes
        class crossSectionMeshGenerator
        {
            // Private data
                //- reference to Time
                const Time& runTime_;

                //- region name
                const fileName regionName_;

                //- time step
                const fileName timeStep_;

                //- name of contact patches
                HashSet<const word> contactPatchNames_;

                //- const pointer to a surface
                const triSurf* surfacePtr_;

                //- const pointer to a modified surface due to anisotropy
                const triSurf* modSurfacePtr_;

                //- const reference to meshDict
                IOdictionary meshDict_;

                //- dictionary user-specific settings
                dictionary userSettingsDict_;

                //- geometry deviation tolerance
                scalar geometryTolerance_;

                //- minimum number of faces in patch
                label minNumFacesInPatch_;

                //- select additional layer of cells near the ones selected
                //- for refinement
                label nAdditionalLayers_;

                //- pointer to the mesh
                polyMeshGen* meshPtr_;

                //- pointer to the octree used for extraction of cartesian mesh
                meshOctree* meshOctreePtr_;

                //- pointer to meshOctree
                meshOctree* octreePtr_;

                //- desired cell size for the contact patch
                direction contactRefLevel_;

                //- minimum allowed cell size for the contact patch
                direction maxContactRefLevel_;

                //- uniform refinement at the boundary
                bool enforceUniformRefinement_;

                //- enforce contact cell size
                bool enforceContactCellSize_;

            // Private member functions
                //- cleanup of small surface patches
                void cleanupSmallPatches();

                //- create the octree from the current input geometry
                void createOctree();

                //- parse user dictionary
                void parseUserDictionary();

                //- create cartesian mesh
                void createCartesianMesh();

                //- prepare mesh surface
                bool surfacePreparation();

                //- map mesh to the surface and untangle surface
                void mapMeshToSurface();

                //- capture edges and corners
                void extractPatches();

                //- capture edges and corners
                void mapEdgesAndCorners();

                //- optimise surface mesh
                void optimiseMeshSurface();

                //- add boundary layers
                void generateBoundaryLayers();

                //- refine boundary layers
                void refBoundaryLayers();

                //- replace boundaries
                void replaceBoundaries();

                //- check if the mesh surface is within tolerance
                //- required by the user
                bool checkGeometryDeviation();

                //- renumber mesh
                void renumberMesh();

                //- execute the 2D workflow
                void generate2DMesh();

                //- disallow bitwise assignment
                void operator=(const crossSectionMeshGenerator&);

                //- disallow copy construct
                crossSectionMeshGenerator
                (
                    const crossSectionMeshGenerator&
                );

        public:

            // Constructors
                //- construct from surface mesh and dictionary
                crossSectionMeshGenerator
                (
                    const triSurf* surfPtr,
                    const IOdictionary& meshDict,
                    const dictionary& userDict,
                    const Time& runTime,
                    const fileName regionName,
                    const fileName timeStep,
                    const HashSet<const word>& contactPatchNames,
                    const bool enforceUniformRefinement = false
                );

            // Destructor

                ~crossSectionMeshGenerator();

            // Public member functions
                //- provide a 2D cross-section mesh
                polyMeshGen* meshPtr() const;
        };

    // Private member functions

        //- Disallow default bitwise copy construct
        rollingMillMesh(const rollingMillMesh&);

        //- Disallow default bitwise assignment
        void operator=(const rollingMillMesh&);

public:

    // Constructors

        //- Construct from time and the region name
        rollingMillMesh
        (
            const Time&,
            const fileName regionName = "",
            const fileName timeStep = ""
        );

    // Destructor

        ~rollingMillMesh();


    // Member Functions
        //- generate a wire mesh
        void generateWire();

        //- generate a mesh for a roller
        void generateRollers();

        //- generate a mesh of a die
        void generateDie();

        //- write the mesh
        void writeMesh() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
