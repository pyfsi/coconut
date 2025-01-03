#include <SMACseApi.h>          // SMACseApi.h is the main include for the CSE API
#include <SMACseApiService.h>
#include <SMACseApiEnums.h>     // SMACseApiEnums.h is included by SMACseApi.h. It is listed for demo purpose only.
#include <SMACseApiMesh.h>
#include <iostream>             // For input/output operations
#include <vector>               // For std::vector
#include <map>                  // For std::map
#include <thread>               // For std::this_thread::sleep_for
#include <chrono>               // For time utilities
#include <fstream>              // For file stream operations
#include <sstream>              // For string stream utilities
#include <algorithm>            // For std::find and other algorithms
#include <cstring>              // For std::strcmp


using namespace std;

// Function arguments prescribed by SCE
void putField(const char* fieldType, const char* meshName, const char* collection, unsigned int nMembers, int* members, void* targetBuffer, void* userData);
void myMessageHandler(const SMACseMsgHandlerSeverity severity, const char* msg, const int numNodes, const int* nodes, const int numElements, const int* elements);

// Own functions
void readInput(const string filenameInput, double& dt, string& port, unsigned int& nDim, unsigned int& nModelParts);
void connectToCSE(void* userData, const string& port);
int findIndex(const vector<unsigned int>& vec, unsigned int target);
vector<double> getElemCentroidCoordinates(const vector<unsigned int>& nodeLabels, const vector<double>& nodeCoordinates, const vector<unsigned int>& connectivity, const unsigned int& nElems, const unsigned int& nodePerElem, const unsigned int& nDim);
int writeRemoteMeshData(const char* localMeshName, unsigned int* nNodes, unsigned int* nElems, const unsigned int& modelPartID);
void readData(const string& filename, vector<double>& values);
template <typename T>
void writeData(const string& filename, const vector<T>& values, const unsigned int& dim, const vector<unsigned int>& labels = {});
bool checkMessage(const string& message);
void waitMessage(const string& message);
void sendMessage(const string& message);

struct MeshInputData {
    string meshName;
    vector<double>* pressure;
    vector<double>* traction;
    unsigned int modelPartID;

    MeshInputData() : pressure(nullptr), traction(nullptr) {}
};

unsigned int timeStep = 0; // Time step counter
unsigned int iteration = 0; // Coupling iteration counter
unsigned int debug = 0; // Debug boolean
unsigned int nDim = 0; // Number of dimensions
const string extension = ".txt"; // Data file extension


int ABQmain(int argc, char** argv){ // Instead of: int main(){ // To be able to compile with abaqus make
    double startTime = 0; // Start time of time step
    double endTime = 0; // End time of time step
    double CSETargetTime; // End time of time step determined by CSE
    double CSEdt; // Time step size determined by CSE
    int lockstep; // Bool indicating requirement to used dt provided by CSE

    double dt; // Time step size
    string port; // Port for connection
    unsigned int nModelParts; // Number of model parts
    const string filenameInput = "AbaqusWrapper_input.txt";
    readInput(filenameInput, dt, port, nDim, nModelParts); // Read input from CoCoNuT wrapper

    const string filenamePressure = "../pressure_mp"; // To be appended with <modelPartID>.txt
    const string filenameTraction = "../traction_mp"; // To be appended with <modelPartID>.txt
    const string filenameDisplacement = "../displacement_mp"; // To be appended with <modelPartID>.txt

    SMACseSolutionStatus status; // Generic status for CSE
    SMACseGetFieldStatus gStatus; // GetField status

    map<string, MeshInputData> meshInputDataMap; // Create a nested map with the structures containing vectors pressure and traction for each model part
    void* const userData = &meshInputDataMap; // Opaque pointer to meshInputDataMap, used to store pressure and traction vectors for each model part accessed in putField

    SMACseInitialize(); // Initialize CSE
    connectToCSE(userData, port); // Establish connection with CSE

    SMACseCreateFieldDefinition("pressure", SMACseFieldPos_atElemCentroid, SMACseFieldAlgType_Scalar, SMACseFieldDataType_Double); // Define field pressure
    SMACseFieldAlgebraicType vectorType;
    switch (nDim) {
        case 2:
            vectorType = SMACseFieldAlgType_Vector2d;
            break;
        case 3:
            vectorType = SMACseFieldAlgType_Vector3d;
            break;
        default:
            cerr << "Unknown number of dimensions " << nDim << endl;
            exit(1);
    }
    SMACseCreateFieldDefinition("traction_vector", SMACseFieldPos_atElemCentroid, vectorType, SMACseFieldDataType_Double); // Define field traction
    SMACseCreateFieldDefinition("displacement", SMACseFieldPos_atNode, vectorType, SMACseFieldDataType_Double); // Define field displacement

    cout << "Fields created" << endl;

    vector<string> meshNameArray(nModelParts); // Mesh name for each model part
    vector<unsigned int> nNodesArray(nModelParts); // Number of nodes for each model part
    vector<unsigned int> nElemsArray(nModelParts); // Number of elements for each model part

    for (unsigned int i = 0; i < nModelParts; ++i) { // Loop over model parts
        meshNameArray[i] = "wrapperMesh" + to_string(i);
        const char* meshName = meshNameArray[i].c_str();
        cout << "Creating mesh and registering fields mesh of model part " << i << " corresponding to " << meshName << endl;

        SMACseMesh* mesh = SMACseCreateMesh(meshName); // Create a NULL mesh

        // Register incoming and outgoing fields
        SMACseRegisterOutgoingField("pressure", meshName); // To Abaqus
        SMACseRegisterOutgoingField("traction_vector", meshName); // To Abaqus
        SMACseRegisterIncomingField("displacement", meshName); // From Abaqus
    }

    SMACseRegisterModel(); // Register model with the NULL meshes and fields with the CSE
    cout << endl << "Model registered" << endl;

    for (unsigned int i = 0; i < nModelParts; ++i) { // Loop over model parts
        meshNameArray[i] = "wrapperMesh" + to_string(i);
        const char* meshName = meshNameArray[i].c_str();
        cout << "> Getting mesh of model part " << i << " corresponding to " << meshName << endl;

        int result = writeRemoteMeshData(meshName, &nNodesArray[i], &nElemsArray[i], i); // Write mesh data and obtain nNodes and nElems
        if (result == 0) {
            cout << "Mesh written successfully" << endl;
        } else {
            cerr << "Something went wrong in getting and writing the remote mesh data" << endl;
            return -1;
        }
    }

    sendMessage("export_mesh_data_ready"); // Notify CoCoNuT mesh has been exported

    vector<vector<double>> pressureArray(nModelParts); // Pressure vector for each model part
    vector<vector<double>> tractionArray(nModelParts); // Traction vector for each model part
    vector<vector<double>> displacementArray(nModelParts); // Displacement vector for each model part

    for (unsigned int i = 0; i < nModelParts; ++i) { // Loop over model parts
        pressureArray[i] = vector<double>(nElemsArray[i]);
        tractionArray[i] = vector<double>(nDim * nElemsArray[i]);
        displacementArray[i] = vector<double>(nDim * nNodesArray[i]);

        MeshInputData meshInputData;
        meshInputData.meshName = meshNameArray[i];
        meshInputData.pressure = &pressureArray[i];
        meshInputData.traction = &tractionArray[i];
        meshInputData.modelPartID = i;
        meshInputDataMap[meshNameArray[i]] = meshInputData;
    }

    status = SMACseInitialConditionsStart(); // Notify CSE of readiness to start initial condition computation
    status = SMACseInitialConditionsEnd(); // Notify CSE that computation of initial conditions ended

    while (true) { // FSI loop
        if (checkMessage("next")) { // Start next step
            startTime = endTime;
            cout << "\n>> Timestep " << ++timeStep << endl;
            endTime = timeStep * dt;
            iteration = 0;

            status = SMACseGetTargetTime(startTime, dt, &CSETargetTime, &CSEdt, &lockstep, NULL); // Gets the preferred time step size from CSE (here constant dt)

            if (abs(endTime - CSETargetTime) > dt * 1e-8) {
                cerr << "The time step size prescribed by the CoSimulation Engine ("
                     << CSEdt << ") doesn't match the set dt (" << dt << ")\n"
                     << "CSETargetTime = " << CSETargetTime << " endTime = " << endTime << endl;
                exit(1);
            }

            status = SMACseNotifyStart(); // Notify readiness to start new time step to CSE
            sendMessage("next_ready");
        }

        if (checkMessage("continue")) { // Start next coupling iteration
            cout << "> Coupling iteration " << ++iteration << endl;

            // Read data from CoCoNuT
            for (unsigned int i = 0; i < nModelParts; ++i) { // Loop over model parts
                readData(filenamePressure + to_string(i) + extension, pressureArray[i]); // Read in pressure
                readData(filenameTraction + to_string(i) + extension, tractionArray[i]); // Read in traction vector
                if (debug) {
                    string extraInfo = + "_ts" + to_string(timeStep) + "_it" + to_string(iteration);
                    writeData(string("read_pressure_mp") + to_string(i) + extraInfo + extension, pressureArray[i], 1); // Write read pressure
                    writeData(string("read_traction_mp") + to_string(i) + extraInfo  + extension, tractionArray[i], nDim); // Write read traction
                }
            }

            // Do iterations
            status = SMACseNotifyIteration(endTime, 0); // Run Abaqus

            // Write data from Abaqus
            for (unsigned int i = 0; i < nModelParts; ++i) { // Loop over model parts
                gStatus = SMACseGetField("displacement", meshNameArray[i].c_str(), NULL, nNodesArray[i], NULL, endTime, static_cast<void*>(displacementArray[i].data()), displacementArray[i].size() * sizeof(double), NULL); // Get displacement from CSE
                // Check gstatus
                if (gStatus == SMACseGetFieldStatus_OrphanedMembers) {
                    cerr << "GetFieldStatus: One or more orphan nodes / elements detected." << endl;
                } else if (gStatus == SMACseGetFieldStatus_SystemError) {
                    cerr << "GetFieldStatus: System level error detected, with details provided to message handler." << endl;
                }
                writeData(filenameDisplacement + to_string(i) + extension, displacementArray[i], nDim); // Write displacement
            }

            sendMessage("continue_ready");
        }

        if (checkMessage("end_step")) { // Stop simulation
            cout << "  Coupling iteration " << iteration << " end" << endl;
            status = SMACseNotifyIteration(endTime, 1); // Stop iterating
            status = SMACseNotifyEnd(endTime); // Notify completion of time step to CSE
            sendMessage("end_step_ready");
        }

        if (checkMessage("save")) { // Save data
            sendMessage("save_ready");
        }

        if (checkMessage("stop")) { // Stop simulation
            cout << "\n>> Stop " << endl;
            sendMessage("stop_ready");
            SMACseDisconnect(); // Terminate connection
            SMACseFinalize(); // Finalize
            break;
        }
    }

    return 0;
}

void putField(const char* fieldType, const char* meshName, const char* collection, unsigned int nMembers, int* members, void* targetBuffer, void* userData) {
    // The putField function fills the provided buffer (by CSE) with results corresponding to the named field

    double* target = static_cast<double*>(targetBuffer); // Cast buffer to double pointer
    map<string, MeshInputData> meshInputDataMap = *static_cast<map<string, MeshInputData>*>(userData); // Cast pointer to original type and dereference

    MeshInputData meshInputData = meshInputDataMap[string(meshName)];
    if (strcmp(fieldType, "pressure") == 0) {
        vector<double> pressure = *meshInputData.pressure;
        copy(pressure.data(), pressure.data() + nMembers, target);
    } else if (strcmp(fieldType, "traction_vector") == 0) {
        vector<double> traction = *meshInputData.traction;
        copy(traction.data(), traction.data() + nDim * nMembers, target);
    }

    return;
}

void myMessageHandler(const SMACseMsgHandlerSeverity severity, const char* msg, const int numNodes, const int* nodes, const int numElements, const int* elements) {
    switch (severity) {

    case SMACseMsgHandler_InformationalMessage:  // -- Informational Message, continue execution
        fprintf(stderr, "%s", msg);
        break;

    case SMACseMsgHandler_WarningMessage:        // -- Warning Message, continue execution
        fprintf(stderr, "%s", msg);
        break;

    case SMACseMsgHandler_ErrorMessage:          // -- Error Message, terminate execution
        fprintf(stderr, "%s", msg);
        exit(1);
        break;
    }
}

void readInput(const string filenameInput, double& dt, string& port, unsigned int& nDim, unsigned int& nModelParts) {
    ifstream file(filenameInput);

    if (!file) {
        cerr << "File " << filenameInput << " not found" << endl;
        return;
    }

    string keyword; // To store the keyword
    // Read until all lines are processed
    while (file >> keyword) {
        if (keyword == "dt") {
            file >> dt;
        } else if (keyword == "number_of_model_parts") {
            file >> nModelParts;
        } else if (keyword == "dimensions") {
            file >> nDim;
        } else if (keyword == "port") {
            file >> port;
        } else if (keyword == "debug") {
            file >> debug;
        } else {
            cerr << "Warning: unknown keyword " << keyword << endl;
        }
    }

    // Check if any of the required values were not read
    bool allValuesRead = true;
    if (dt == 0.0) {
        cerr << "Error: 'dt' not read correctly from file" << endl;
        allValuesRead = false;
    }
    if (nModelParts == 0) {
        cerr << "Error: 'number_of_model_parts' not read correctly from file" << endl;
        allValuesRead = false;
    }
    if (nDim == 0) {
        cerr << "Error: 'dimensions' not read correctly from file" << endl;
        allValuesRead = false;
    }
    if (port.empty()) {
        cerr << "Error: 'port' not read correctly from file" << endl;
        allValuesRead = false;
    }
    if (!allValuesRead) {
        cerr << "Error: Some values were not read correctly. Please check the input file." << endl;
        exit(1);
    }

    // Print the values to verify that they were read correctly
    cout << "Input values read from " << filenameInput << endl;;
    cout << "dt: " << dt << endl;
    cout << "nDim: " << nDim << endl;
    cout << "nModelParts: " << nModelParts << endl;
    cout << "port: " << port << endl;
    cout << "debug: " << debug << endl;
}

void connectToCSE(void* userData, const string& port) {
    // Establish connection between CSE Director and all other clients participating in the co-simulation.

    const char* codeName = "AbaqusWrapper"; // Must match the code name for component in CSE configuration file
    const char* jobName = "AbaqusWrapper"; // Must match component instance name in CSE configuration file
    const char* workDir = "./"; // Current directory
    const char* cseDirector = ("localhost:" + port).c_str(); // Port for connection
    const int timeout = 86400; // in seconds
    double startTime = 0.0; // Not yet implemented in CSE

    SMACseConnectStatus cStatus;

    cStatus = SMACseConnect(codeName, jobName, workDir, cseDirector, timeout, startTime, putField, userData, myMessageHandler);

    // Check connection
    if (cStatus == SMACseConnectStatus_failure) {
        myMessageHandler(SMACseMsgHandler_ErrorMessage, "An error occurred during the connection process.\n", 0, NULL, 0, NULL);
        exit(1);
    }
    else {
        printf("Successfully connected with CSE Director and other client process.\n");
    }
}

int findIndex(const vector<unsigned int>& vec, unsigned int target) {
    // Use find to get an iterator to the target
    auto it = find(vec.begin(), vec.end(), target);

    // If the element is found, return the index
    if (it != vec.end()) {
        return distance(vec.begin(), it); // distance gives the index
    } else {
        return -1; // Return -1 if the element is not found
    }
}

vector<double> getElemCentroidCoordinates(const vector<unsigned int>& nodeLabels, const vector<double>& nodeCoordinates, const vector<unsigned int>& connectivity, const unsigned int& nElems, const unsigned int& nodePerElem, const unsigned int& nDim) {
    vector<double> elemCentroidCoordinates(nElems * nDim);
    double xCoordinate;
    double yCoordinate;
    double zCoordinate;
    unsigned int nodeLabel;
    int nodeIndex;

    for (unsigned int i = 0; i < nElems; ++i) { // Loop over elements
        xCoordinate = 0;
        yCoordinate = 0;
        zCoordinate = 0;
        for (unsigned int j = 0; j < nodePerElem; ++j) { // Loop over nodes of element i
            nodeLabel = connectivity[i * nodePerElem + j];
            nodeIndex = findIndex(nodeLabels, nodeLabel);
            if (nodeIndex == -1) {
                cerr << "No node with label " << nodeLabel << endl;
            }
            xCoordinate += nodeCoordinates[nodeIndex * nDim];
            yCoordinate += nodeCoordinates[nodeIndex * nDim + 1];
            zCoordinate += nodeCoordinates[nodeIndex * nDim + 2];
        }
        elemCentroidCoordinates[i * nDim] = xCoordinate / nodePerElem;
        elemCentroidCoordinates[i * nDim + 1] = yCoordinate / nodePerElem;
        elemCentroidCoordinates[i * nDim + 2] = zCoordinate / nodePerElem;
    }

    return elemCentroidCoordinates;
}

int writeRemoteMeshData(const char* localMeshName, unsigned int* nNodesPtr, unsigned int* nElemsPtr, const unsigned int& modelPartID) {
    cout << "Retrieving remote mesh data for mesh " << localMeshName << endl;

    unsigned int nNodes;
    unsigned int nElems;

    unsigned int nElemCollections;
    unsigned int nDimLocal;
    SMACseMesh* remoteMesh = SMACseGetRemoteMesh(localMeshName, &nNodes, &nElemCollections, &nDimLocal); // Get remote mesh (CSE mapper is bypassed)

    if (nDimLocal != nDim) {
        cerr << "Error: Number of dimensions doesn't match the number of dimensions read from input" << endl;
        return 1;
    }

    if (nElemCollections == 0) {
        cerr << "Error: No element collections in Abaqus mesh, verify that the co-simulation interface region is a surface" << endl;
        return 1;
    }
    else if (nElemCollections > 1) {
        cerr << "Error: More than one element collection in Abaqus mesh, this is currently not supported" << endl;
        return 1;
    }

    vector<unsigned int> nodeLabels(nNodes);
    vector<double> nodeCoordinates(nDim * nNodes);
    SMACseGetNodeCollection(remoteMesh, (unsigned int*)nodeLabels.data(), (double*)nodeCoordinates.data());

    unsigned int nodePerElem;
    const char* elementType = SMACseQueryElemCollection(remoteMesh, 0, &nElems, &nodePerElem);

    if (!strstr(elementType, "SMACseSurfTopo")) {
        cerr << "Error: Abaqus elements are not surface elements" << endl;
        return 1;
    }

    vector<unsigned int> elemLabels(nElems);
    vector<unsigned int> connectivity(nElems * nodePerElem);
    SMACseGetElemCollection(remoteMesh, 0, (unsigned int*)elemLabels.data(), (unsigned int*)connectivity.data());

    vector<double> elemCentroidCoordinates = getElemCentroidCoordinates(nodeLabels, nodeCoordinates, connectivity, nElems, nodePerElem, nDim);

    cout << "    number of dimensions " << nDimLocal << endl;
    cout << "    number of nodes " << nNodes << endl;
    cout << "    number of elements " << nElems << endl;
    cout << "    number of element collections " << nElemCollections << " (1 is only allowed value, handling more than one element collection is not yet implemented)"<< endl;
    cout << "    element type " << elementType << endl;
    cout << "    nodes per element " << nodePerElem << endl;

    if (debug) {
        writeData("../node_labels_mp" + to_string(modelPartID) + extension, nodeLabels, 1);
        writeData("../element_labels_mp" + to_string(modelPartID) + extension, elemLabels, 1);
        writeData("../connectivity_mp" + to_string(modelPartID) + extension, connectivity, nodePerElem);
    }
    writeData("../initial_node_coordinates_mp" + to_string(modelPartID) + extension, nodeCoordinates, nDim, nodeLabels);
    writeData("../initial_element_centroid_coordinates_mp" + to_string(modelPartID) + extension, elemCentroidCoordinates, nDim, elemLabels);

    *nNodesPtr = nNodes;
    *nElemsPtr = nElems;

    return 0;
}

string stripFilename(const string& path) {
    // Find the position of the last slash (to remove the path)
    size_t lastSlashPos = path.find_last_of("/\\");
    std::string filename = (lastSlashPos == std::string::npos) ? path : path.substr(lastSlashPos + 1);

    // Find the position of the last dot (to remove the extension)
    size_t lastDotPos = filename.find_last_of(".");
    std::string baseFilename = (lastDotPos == std::string::npos) ? filename : filename.substr(0, lastDotPos);

    return baseFilename;
//    // In c++17 ( #include <filesystem> )
//    filesystem::path p(filename); // Convert the string into a std::filesystem::path object
//    return p.stem().string(); // Get the filename without directory without extension
}

void readData(const string& filename, vector<double>& values) {
    ifstream file(filename);  // Open the file for reading

    if (!file) {
        cerr << "Error: Could not open file " << filename << endl;
        return;
    }

    size_t index = 0; // Keep track of the current position in the vector
    string line;

    // Skip header (2 lines)
    getline(file, line);
    getline(file, line);

    while (getline(file, line)) { // Read the file line by line
        istringstream iss(line);  // Parse the line
        double num;

        while (iss >> num) { // Extract float numbers
            if (index >= values.size()) {
                cerr << "Error: File contains more data than allocated in the vector." << endl;
                return;
            }
            values[index++] = num; // Fill the vector at the current index
        }
    }

    if (index < values.size()) {
        cerr << "Warning: Not all elements in the vector were filled." << endl;
    }

    file.close();
}

template <typename T>
void writeData(const string& filename, const vector<T>& values, const unsigned int& dim, const vector<unsigned int>& labels) {
    // filename is the name of the file to be written to
    // values is the vector that should be written, passed as reference
    // labels is a vector containing optional labels (use NULL to disable)
    // dim is the dimension (number of columns to be written)

    ofstream file(filename); // Open the file for writing

    if (!file) {
        cerr << "Error: Could not open file " << filename << " for writing." << endl;
        return;
    }

    // Print header
    file << "# " << stripFilename(filename) << ": timestep " << timeStep << ", iteration " << iteration << " #" << values.size() / dim << endl;
    string labelName;
    if (!labels.empty()) {
        labelName = "label";
    }
    string valName;
    switch (dim) {
        case 2:
            valName = "x y";
            break;
        case 3:
            valName = "x y z";
            break;
        default:
            break;
    }
    file << "# " << labelName << " " << valName << endl;

    for (size_t i = 0; i < values.size(); ++i) {
        // Print label on each row
        if (!labels.empty() && (i % dim) == 0) {
            file << labels[i / dim] << " ";
        }

        // Print value
        file << values[i];

        // Add a space after each number, and a newline after every 'dimension' numbers
        if ((i + 1) % dim == 0) {
            file << "\n";
        } else {
            file << " ";
        }
    }

    file.close(); // Close the file after writing
}

bool checkMessage(const string& message) {
    string filename = "../" + message + ".coco";
    ifstream file(filename);
    if (file.good()) { // Returns true if the file is accessible
        if (debug) {
            cout << "Received coconut message " << message << endl;
        }
        remove(filename.c_str());
        return file.good();
    } else {
        return false;
    }
}

void waitMessage(const string& message) {
    while (true) {
        if (checkMessage(message)) {
            return;
        }
        this_thread::sleep_for(chrono::milliseconds(100));
    }
}

void sendMessage(const string& message) {
    ofstream file("../" + message + ".coco");
    if (!file) {
        cerr << "Error creating file: " << message << endl;
        return;
    }
    if (debug) {
        cout << "Sent coconut message " << message << endl;
    }
    file.close();
}
