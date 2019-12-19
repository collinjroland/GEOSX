/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior
 * University Copyright (c) 2018-2019 Total, S.A Copyright (c) 2019-     GEOSX
 * Contributors All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS
 * files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file VTKWriter.cpp
 */

#include "VTKWriter.hpp"
#include <sys/stat.h>

#include "dataRepository/Wrapper.hpp"
#include "managers/DomainPartition.hpp"

#ifdef GEOSX_USE_VTK
#include "vtkPoints.h"
#include "vtkPolyhedron.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXMLPUnstructuredGridWriter.h"
#endif

namespace geosx {
using namespace dataRepository;

namespace {
/// Map from GEOSX type to VTK cell types
std::unordered_map<string, string> geosxToVTKCellTypeMap = {
    {"C3D4", "VTK_TETRA"},      // 4 nodes
    {"C3D5", "VTK_PYRAMID"},    // 5 nodes
    {"C3D6", "VTK_WEDGE"},      // 6 nodes
    {"C3D8", "VTK_HEXAHEDRON"}, // 8 nodes
    {"", "VTK_QUAD"},           // 4 nodes, 2D
    {"POLY", "VTK_POLYHEDRON"}  // undefined number of nodes
};

/// Map from GEOSX type to VTK cell types
std::unordered_map<string, int> geosxToVTKCellTypeIdMap = {
    {"C3D4", 10}, // tetra
    {"C3D5", 14}, // pyramid
    {"C3D6", 13}, // wedge
    {"C3D8", 12}, // hexahedron
    {"", 9},      // quad
    {"POLY", 42}  // polyhedron
};

std::unordered_map<std::type_index, string> geosxToVTKTypeMap = {
    {std::type_index(typeid(integer)), "Int32"},
    {std::type_index(typeid(localIndex)), "Int64"},
    {std::type_index(typeid(globalIndex)), "Int64"},
    {std::type_index(typeid(real32)), "Float32"},
    {std::type_index(typeid(real64)), "Float64"},
    {std::type_index(typeid(r1_array)), "Float64"},
    {std::type_index(typeid(real64_array)), "Float64"},
    {std::type_index(typeid(real64_array2d)), "Float64"},
    {std::type_index(typeid(real64_array3d)), "Float64"},
    {std::type_index(typeid(real32_array)), "Float32"},
    {std::type_index(typeid(real32_array2d)), "Float32"},
    {std::type_index(typeid(real32_array3d)), "Float32"},
    {std::type_index(typeid(integer_array)), "Int32"},
    {std::type_index(typeid(localIndex_array)), "Int64"},
    {std::type_index(typeid(localIndex_array2d)), "Int64"},
    {std::type_index(typeid(localIndex_array3d)), "Int64"},
    {std::type_index(typeid(globalIndex_array)), "Int64"},
    {std::type_index(typeid(globalIndex_array2d)), "Int64"},
    {std::type_index(typeid(globalIndex_array3d)), "Int64"}};
} /* namespace */

vtkSmartPointer<vtkPoints> extractPoints(NodeManager const *nodeManager)
{
  int n_vertices = nodeManager->size(); // check type output
  r1_array const &vertices = nodeManager->referencePosition();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(
      n_vertices); // pre-allocate memory and set maximum Id
  vtkIdType id[1];
  for (int i = 0; i < n_vertices; ++i) {
    id[0] = i;
    auto vertex = vertices[i];
    // Insert vertex in vtkPoints array
    points->SetPoint(id, vertex[0], vertex[1], vertex[2]);
  }

  points->Squeeze(); // Probably unnecessary, since we're allocating exactly the
                     // right number of points
  return points
}

vtkSmartPointer<vtkCellArray>
extractConnectivity(ElementRegionManager const *elemManager, int *cellTypes)
{
  localIndex nCells = elemManager->getNumberOfElements<CellElementSubRegion>();
  int currentIndex = 0;

  vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
  cells->SetNumberOfCells(nCells); // May need nCells as vtkIdType

  elemManager->forElementRegionsComplete<CellElementRegion>(
      [&](localIndex const GEOSX_UNUSED_ARG(er), auto const *const elemRegion) {
        elemRegion->template forElementSubRegions<CellElementSubRegion>(
            [&](auto const *const elemSubRegion) {
              int type = geosxToVTKCellTypeIdMap.at(
                  elemSubRegion->GetElementTypeString());
              for (localIndex k = 0; k < elemSubRegion->size(); ++k) {
                cellTypes[currentIndex] = type;
                currentIndex++;
                auto &connectivities = elemSubRegion->nodeList(k);
                // cases for each type
                if (type == 12) {
                  // special loop for hexahedra
                  vtkSmartPointer<vtkHexahedron> currentCell =
                      vtkSmartPointer<vtkHexahedron>::New();
                  currentCell->GetPointIds->SetId(0, connectivities[0]);
                  currentCell->GetPointIds->SetId(1, connectivities[1]);
                  currentCell->GetPointIds->SetId(2, connectivities[3]);
                  currentCell->GetPointIds->SetId(3, connectivities[2]);
                  currentCell->GetPointIds->SetId(4, connectivities[4]);
                  currentCell->GetPointIds->SetId(5, connectivities[5]);
                  currentCell->GetPointIds->SetId(6, connectivities[7]);
                  currentCell->GetPointIds->SetId(7, connectivities[6]);
                  cells->InsertNextCell(currentCell);
                } else {
                  cells->InsertNextCell(connectivities.size());
                  for (int i = 0; i < connectivities.size(); ++i) {
                    // May need conversion to vtkIdType
                    cells->InsertCellPoint(connectivities[i]);
                    // Probably needs third option with face stream for
                    // arbitrary polyhedra
                  }
                }
              }
            });
      });

  cells->Squeeze();
  return cells;
}

void VTKWriter::Write(double timeStep, DomainPartition const &domain)
{
  ElementRegionManager const *elemManager =
      domain.getMeshBody(0)->getMeshLevel(0)->getElemManager();
  NodeManager const *nodeManager =
      domain.getMeshBody(0)->getMeshLevel(0)->getNodeManager();
  auto dataView =
      elemManager->ConstructViewAccessor<cType>(std::get<0>(cellField));
  localIndex totalNumberOfCells =
      elemManager->getNumberOfElements<CellElementSubRegion>();

  // Extract data from mesh
  vtkSmartPointer<vtkPoints> points = extractPoints(nodeManager);
  int cellTypes[totalNumberOfCells];
  vtkSmartPointer<vtkCellArray> cellArray =
      extractConnectivity(elemManager, cellTypes);
  // split into Point and Cell data fields?
  auto fields = extractCellData(dataView, elemManager);

  // Put data into VTK unstructured grid
  vtkSmartPointer<vtkUnstructuredGrid> ugrid =
      vtkSmartPointer<vtkUnstructuredGrid>::New();
  ugrid->SetPoints(points);
  ugrid->SetCells(cellTypes, cellArray);
  ugrid->; // set field data

  // Write results
  vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer =
      vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
  writer->SetInputData(ugrid);
  writer->SetFileName(outFile.c_str());
  writer->SetNumberOfPieces(nPieces);
  writer->SetStartPiece(0);
  writer->SetEndPiece(nPieces - 1);
  if (!binary):
    writer->SetDataModeToAscii();
  writer.Update();
}



/////////////////////////////////////////////////////
// FUNCTIONS BELOW THIS COMMENT FOR REFERENCE ONLY //
/////////////////////////////////////////////////////



template <typename T>
void WriteCellAsciiData(
    ElementRegionManager::ElementViewAccessor<T> const &dataView,
    ElementRegionManager const *const elemManager)
{
  elemManager->forElementRegionsComplete<CellElementRegion>(
      [&](localIndex const er, auto const *const elemRegion) {
        elemRegion->template forElementSubRegionsIndex<CellElementSubRegion>(
            [&](localIndex const esr, auto const *const elemSubRegion) {
              for (localIndex ei = 0; ei < elemSubRegion->size(); ei++) {
                m_outFile << dataView[er][esr][ei] << "\n";
              }
            });
      });
}

template <typename T> void WriteNodeAsciiData(Wrapper<T> const &dataView)
{
  auto &viewRef = dataView.reference();
  for (localIndex i = 0; i < viewRef.size(); i++) {
    m_outFile << viewRef[i] << "\n";
  }
}

class CustomVTUXMLWriter {
public:
  /*!
   * @brief Write the vertices coordinates
   * @param[in] vertices table of vertice coordinates
   * @param[in] binary tells wether or not the data should be written in binary
   * format
   */
  void WriteVertices(r1_array const &vertices, bool binary)
  {
    if (binary) {
      WriteBinaryVertices(vertices);
    } else {
      WriteAsciiVertices(vertices);
    }
  }

  /*!
   * @brief Write the cell connectivities
   */
  void WriteCellConnectivities(ElementRegionManager const *const elemManager,
                               bool binary)
  {
    if (binary) {
      localIndex totalNumberOfConnectivities = 0;
      elemManager->forElementRegionsComplete<CellElementRegion>(
          [&](localIndex const GEOSX_UNUSED_ARG(er),
              auto const *const elemRegion) {
            elemRegion->template forElementSubRegions<CellElementSubRegion>(
                [&](auto const *const elemSubRegion) {
                  totalNumberOfConnectivities +=
                      elemSubRegion->size() *
                      elemSubRegion->numNodesPerElement();
                });
          });
      WriteSize(totalNumberOfConnectivities, sizeof(localIndex));
      WriteBinaryConnectivities(elemManager);
    } else {
      WriteAsciiConnectivities(elemManager);
    }
  }

  /*!
   * @brief Write the offsets
   * @details for a full hex mesh : 0, 8, 16, 24....
   */
  void WriteCellOffsets(ElementRegionManager const *const elemManager,
                        bool binary)
  {
    if (binary) {
      WriteSize(elemManager->getNumberOfElements<CellElementSubRegion>(),
                sizeof(localIndex));
      WriteBinaryOffsets(elemManager);
    } else {
      WriteAsciiOffsets(elemManager);
    }
  }

  /*!
   * @brief Write the Cell types
   */
  void WriteCellTypes(ElementRegionManager const *const elemManager,
                      bool binary)
  {
    if (binary) {
      WriteSize(elemManager->getNumberOfElements<CellElementSubRegion>(),
                sizeof(integer));
      WriteBinaryTypes(elemManager);
    } else {
      WriteAsciiTypes(elemManager);
    }
  }

  template <typename T>
  void
  WriteCellData(ElementRegionManager::ElementViewAccessor<T> const &dataView,
                ElementRegionManager const *const elemManager, bool binary)
  {
    if (binary) {
      WriteCellBinaryData(dataView, elemManager);
    } else {
      WriteCellAsciiData(dataView, elemManager);
    }
  }

  template <typename T>
  void WriteNodeData(Wrapper<T> const &dataView, bool binary)
  {
    if (binary) {
      WriteNodeBinaryData(dataView);
    } else {
      WriteNodeAsciiData(dataView);
    }
  }

private:
  void WriteAsciiVertices(r1_array const &vertices)
  {
    for (auto vertex : vertices) {
      m_outFile << vertex << "\n";
    }
  }

  void WriteAsciiConnectivities(ElementRegionManager const *const elemManager)
  {
    elemManager->forElementRegionsComplete<CellElementRegion>(
        [&](localIndex const GEOSX_UNUSED_ARG(er),
            auto const *const elemRegion) {
          elemRegion->template forElementSubRegions<CellElementSubRegion>(
              [&](auto const *const elemSubRegion) {
                integer type = geosxToVTKCellTypeMap.at(
                    elemSubRegion->GetElementTypeString());
                auto &connectivities = elemSubRegion->nodeList();
                if (type == 12) // Special case for hexahedron because of the
                                // internal ordering
                {
                  for (localIndex i = 0; i < connectivities.size() / 8; i++) {
                    m_outFile << connectivities[i][0] << " ";
                    m_outFile << connectivities[i][1] << " ";
                    m_outFile << connectivities[i][3] << " ";
                    m_outFile << connectivities[i][2] << " ";
                    m_outFile << connectivities[i][4] << " ";
                    m_outFile << connectivities[i][5] << " ";
                    m_outFile << connectivities[i][7] << " ";
                    m_outFile << connectivities[i][6] << " ";
                    m_outFile << "\n";
                  }
                } else {
                  for (localIndex i = 0; i < connectivities.size(); i++) {
                    m_outFile << connectivities.data()[i] << " ";
                  }
                  m_outFile << "\n";
                }
              });
        });
  }

  void WriteAsciiOffsets(ElementRegionManager const *const elemManager)
  {
    localIndex curOffset =
        elemManager->GetRegion(0)->GetSubRegion(0)->numNodesPerElement();
    elemManager->forElementRegionsComplete<CellElementRegion>(
        [&](localIndex const GEOSX_UNUSED_ARG(er),
            auto const *const elemRegion) {
          elemRegion->template forElementSubRegions<CellElementSubRegion>(
              [&](auto const *const elemSubRegion) {
                localIndex offSetForOneCell =
                    elemSubRegion->numNodesPerElement();
                for (localIndex i = 0; i < elemSubRegion->size(); i++) {
                  m_outFile << curOffset << "\n";
                  curOffset += offSetForOneCell;
                }
              });
        });
  }

  void WriteAsciiTypes(ElementRegionManager const *const elemManager)
  {
    elemManager->forElementRegionsComplete<CellElementRegion>(
        [&](localIndex const GEOSX_UNUSED_ARG(er),
            auto const *const elemRegion) {
          elemRegion->template forElementSubRegions<CellElementSubRegion>(
              [&](auto const *const elemSubRegion) {
                integer type = geosxToVTKCellTypeMap.at(
                    elemSubRegion->GetElementTypeString());
                for (localIndex i = 0; i < elemSubRegion->size(); i++) {
                  m_outFile << type << "\n";
                }
              });
        });
  }

  template <typename T>
  void WriteCellBinaryData(
      ElementRegionManager::ElementViewAccessor<T> const &dataView,
      ElementRegionManager const *const elemManager)
  {
    std::stringstream stream;
    WriteSize(elemManager->getNumberOfElements<CellElementSubRegion>(),
              sizeof(dataView[0][0][0]));
    integer multiplier = FindMultiplier(
        sizeof(dataView[0][0][0])); // We do not write all the data at once to
                                    // avoid creating a big table each time.
    string outputString;
    outputString.resize(
        FindBase64StringLength(sizeof(dataView[0][0][0]) * multiplier));
    T dataFragment(multiplier);
    integer countDataFragment = 0;
    elemManager->forElementRegionsComplete<CellElementRegion>(
        [&](localIndex const er, auto const *const elemRegion) {
          elemRegion->template forElementSubRegionsIndex<CellElementSubRegion>(
              [&](localIndex const esr, auto const *const elemSubRegion) {
                if (dataView[er][esr].size() > 0)
                  for (localIndex ei = 0; ei < elemSubRegion->size(); ei++) {
                    dataFragment[countDataFragment++] = dataView[er][esr][ei];
                    if (countDataFragment == multiplier) {
                      stream << stringutilities::EncodeBase64(
                          reinterpret_cast<const unsigned char *>(
                              dataFragment.data()),
                          outputString,
                          sizeof(dataView[0][0][0]) * countDataFragment);
                      countDataFragment = 0;
                    }
                  }
                else {
                  real64_array nanArray(3);
                  nanArray[0] = nanArray[1] = nanArray[2] = std::nan("0");
                  for (localIndex ei = 0; ei < elemSubRegion->size(); ei++) {
                    stream << stringutilities::EncodeBase64(
                        reinterpret_cast<const unsigned char *>(
                            nanArray.data()),
                        outputString, sizeof(real64) * 3);
                  }
                }
              });
        });
    outputString.resize(FindBase64StringLength(sizeof(dataView[0][0][0]) *
                                               (countDataFragment)));
    stream << stringutilities::EncodeBase64(
        reinterpret_cast<const unsigned char *>(dataFragment.data()),
        outputString, sizeof(dataView[0][0][0]) * (countDataFragment));
    DumpBuffer(stream);
  }

  template <typename T> void WriteNodeBinaryData(Wrapper<T> const &dataView)
  {
    std::stringstream stream;
    auto &viewRef = dataView.reference();
    integer multiplier = FindMultiplier(
        sizeof(viewRef[0])); // We do not write all the data at once to avoid
                             // creating a big table each time.
    string outputString;
    outputString.resize(
        FindBase64StringLength(sizeof(viewRef[0]) * multiplier));
    T dataFragment(multiplier);
    integer countDataFragment = 0;
    WriteSize(viewRef.size(), sizeof(viewRef[0]));
    for (localIndex i = 0; i < viewRef.size(); i++) {
      dataFragment[countDataFragment++] = viewRef[i];
      if (countDataFragment == multiplier) {
        stream << stringutilities::EncodeBase64(
            reinterpret_cast<const unsigned char *>(dataFragment.data()),
            outputString, sizeof(viewRef[0]) * countDataFragment);
        countDataFragment = 0;
      }
    }
    outputString.resize(
        FindBase64StringLength(sizeof(viewRef[0]) * (countDataFragment)));
    stream << stringutilities::EncodeBase64(
        reinterpret_cast<const unsigned char *>(dataFragment.data()),
        outputString, sizeof(viewRef[0]) * (countDataFragment));
    DumpBuffer(stream);
  }

  /*!
   * @brief This function is used to compute the minimum number of value of a
   * certain type that can be continously encoded into a base64 to be properly
   * written into the VTU file.
   */
  integer FindMultiplier(integer typeSize)
  {
    integer multiplier = 1;
    while ((multiplier * typeSize) % 6) {
      multiplier++;
    }
    return multiplier;
  }

  integer FindBase64StringLength(integer dataSize)
  {
    integer base64StringLength = (dataSize * 8) / 6;
    while (base64StringLength % 4) {
      base64StringLength++;
    }
    return base64StringLength;
  }

  void DumpBuffer(std::stringstream const &stream)
  {
    m_outFile << stream.rdbuf() << '\n';
  }

private:
  /// vtu output file
  std::ofstream m_outFile;

  /// Space counter to have well indented XML file
  int m_spaceCount;
};

template <>
inline void CustomVTUXMLWriter::WriteCellBinaryData(
    ElementRegionManager::ElementViewAccessor<r1_array> const &dataView,
    ElementRegionManager const *const elemManager)
{
  std::stringstream stream;
  string outputString;
  outputString.resize(FindBase64StringLength(sizeof(real64) * 3));
  WriteSize(elemManager->getNumberOfElements<CellElementSubRegion>() * 3,
            sizeof(real64));
  elemManager->forElementRegionsComplete<CellElementRegion>(
      [&](localIndex const er, auto const *const elemRegion) {
        elemRegion->template forElementSubRegionsIndex<CellElementSubRegion>(
            [&](localIndex const esr, auto const *const elemSubRegion) {
              if (dataView[er][esr].size() > 0)
                for (localIndex ei = 0; ei < elemSubRegion->size(); ei++) {
                  stream << stringutilities::EncodeBase64(
                      reinterpret_cast<const unsigned char *>(
                          dataView[er][esr][ei].Data()),
                      outputString, sizeof(real64) * 3);
                }
              else {
                real64_array nanArray(3);
                nanArray[0] = nanArray[1] = nanArray[2] = std::nan("0");
                for (localIndex ei = 0; ei < elemSubRegion->size(); ei++) {
                  stream << stringutilities::EncodeBase64(
                      reinterpret_cast<const unsigned char *>(nanArray.data()),
                      outputString, sizeof(real64) * 3);
                }
              }
            });
      });
  DumpBuffer(stream);
}

template <>
inline void
CustomVTUXMLWriter::WriteNodeBinaryData(Wrapper<r1_array> const &dataView)
{
  std::stringstream stream;
  auto &viewRef = dataView.reference();
  string outputString;
  outputString.resize(FindBase64StringLength(sizeof(real64) * 3));
  WriteSize(viewRef.size() * 3, sizeof(real64));
  for (localIndex i = 0; i < viewRef.size(); i++) {
    stream << stringutilities::EncodeBase64(
        reinterpret_cast<const unsigned char *>(viewRef[i].Data()),
        outputString, sizeof(real64) * 3);
  }
  DumpBuffer(stream);
}

VTKWriter::VTKWriter(string const &name) : m_baseName(name), m_binary(false)
{
  int const mpiRank = MpiWrapper::Comm_rank(MPI_COMM_GEOSX);
  if (mpiRank == 0) {
    // Declaration of XML version
    auto declarationNode = m_rootFile.append_child(pugi::node_declaration);
    declarationNode.append_attribute("version") = "1.0";

    // Declaration of the node VTKFile
    auto vtkFileNode = m_rootFile.append_child("VTKFile");
    vtkFileNode.append_attribute("type") = "Collection";
    vtkFileNode.append_attribute("version") = "0.1";
    // vtkFileNode.append_attribute("byteOrder") = "LittleEndian";
    // vtkFileNode.append_attribute("compressor") = "vtkZLibDataCompressor";

    // Declaration of the node Collection
    vtkFileNode.append_child("Collection");
    mode_t mode = 0733;
    mkdir(name.c_str(), mode);

    string pvdFileName = name + ".pvd";
    m_rootFile.save_file(pvdFileName.c_str());
  }
}

void VTKFile::Write(double const timeStep, DomainPartition const &domain)
{
  int const mpiRank = MpiWrapper::Comm_rank(MPI_COMM_GEOSX);
  int const mpiSize = MpiWrapper::Comm_size(MPI_COMM_GEOSX);
  ElementRegionManager const *elemManager =
      domain.getMeshBody(0)->getMeshLevel(0)->getElemManager();
  NodeManager const *nodeManager =
      domain.getMeshBody(0)->getMeshLevel(0)->getNodeManager();
  string timeStepFolderName = m_baseName + "/" + std::to_string(timeStep);
  string format;
  if (m_binary) {
    format = "binary";
  } else {
    format = "ascii";
  }

  std::set<std::tuple<string, string, integer, rtTypes::TypeIDs>>
      cellFields; // First : field name, Second : type, Third : field dimension;
  // Find all cell fields to export
  elemManager->forElementRegionsComplete<CellElementRegion>(
      [&](localIndex const GEOSX_UNUSED_ARG(er), auto const *const elemRegion) {
        elemRegion->forElementSubRegions([&](auto const *const subRegion) {
          for (auto const &wrapperIter : subRegion->wrappers()) {
            WrapperBase const *const wrapper = wrapperIter.second;

            if (wrapper->getPlotLevel() < m_plotLevel) {
              // the field name is the key to the map
              string const fieldName = wrapper->getName();
              std::type_info const &typeID = wrapper->get_typeid();
              rtTypes::TypeIDs fieldType =
                  rtTypes::typeID(wrapper->get_typeid());
              if (!geosxToVTKTypeMap.count(typeID)) continue;
              int dimension = 0;
              if (fieldType == rtTypes::TypeIDs::r1_array_id) {
                dimension = 3;
              } else {
                dimension = 1;
              }
              cellFields.insert(std::make_tuple(fieldName,
                                                geosxToVTKTypeMap.at(typeID),
                                                dimension, fieldType));
            }
          }
        });
      });

  std::set<std::tuple<string, string, integer, rtTypes::TypeIDs>>
      nodeFields; // First : field name, Second : type, Third : field dimension;
  // Find all node fields to export
  for (auto const &wrapperIter : nodeManager->wrappers()) {
    WrapperBase const *const wrapper = wrapperIter.second;
    if (wrapper->getPlotLevel() < m_plotLevel) {
      string const fieldName = wrapper->getName();
      std::type_info const &typeID = wrapper->get_typeid();
      if (!geosxToVTKTypeMap.count(typeID)) continue;
      int dimension = 0;
      rtTypes::TypeIDs fieldType = rtTypes::typeID(wrapper->get_typeid());
      if (fieldType == rtTypes::TypeIDs::r1_array_id) {
        dimension = 3;
      } else {
        dimension = 1;
      }
      nodeFields.insert(std::make_tuple(fieldName, geosxToVTKTypeMap.at(typeID),
                                        dimension, fieldType));
    }
  }
  if (mpiRank == 0) {
    /// Add the new entry to the pvd root file
    auto collectionNode = m_rootFile.child("VTKFile").child("Collection");
    auto dataSetNode = collectionNode.append_child("DataSet");
    dataSetNode.append_attribute("timestep") = std::to_string(timeStep).c_str();
    dataSetNode.append_attribute("group") = "";
    dataSetNode.append_attribute("part") = "0";
    string pvtuFileName = timeStepFolderName + "/root.pvtu";
    dataSetNode.append_attribute("file") = pvtuFileName.c_str();

    /// Create the pvtu file for this time step
    // Create a directory for this time step
    mode_t mode = 0733;
    mkdir(timeStepFolderName.c_str(), mode);
    pugi::xml_document pvtuFile;

    // Declaration of XML version
    auto declarationNode = pvtuFile.append_child(pugi::node_declaration);
    declarationNode.append_attribute("version") = "1.0";

    // Declaration of the node VTKFile
    auto vtkFileNode = pvtuFile.append_child("VTKFile");
    vtkFileNode.append_attribute("type") = "PUnstructuredGrid";
    vtkFileNode.append_attribute("version") = "0.1";
    if (m_binary) {
      vtkFileNode.append_attribute("byteOrder") = "LittleEndian";
    }

    // Declaration of the node PUnstructuredGrid
    auto pUnstructureGridNode = vtkFileNode.append_child("PUnstructuredGrid");
    pUnstructureGridNode.append_attribute("GhostLevel") = "1";

    // Declaration the node PPoints
    auto pPointsNode = pUnstructureGridNode.append_child("PPoints");
    // .... and the data array containg the positions
    CreatePDataArray(pPointsNode,
                     geosxToVTKTypeMap.at(std::type_index(typeid(real64))),
                     "Position", 3, format);

    // Declare all the point fields
    auto pPointDataNode = pUnstructureGridNode.append_child("PPointData");
    for (auto &nodeField : nodeFields) {
      CreatePDataArray(pPointDataNode, std::get<1>(nodeField),
                       std::get<0>(nodeField), std::get<2>(nodeField), format);
    }

    // Declaration of the node PCells
    auto pCellsNode = pUnstructureGridNode.append_child("PCells");
    // .... and its data array defining the connectivities, types, and offsets
    CreatePDataArray(pCellsNode,
                     geosxToVTKTypeMap.at(std::type_index(typeid(localIndex))),
                     "connectivity", 1, format); // TODO harcoded for the moment
    CreatePDataArray(pCellsNode,
                     geosxToVTKTypeMap.at(std::type_index(typeid(localIndex))),
                     "offsets", 1, format);
    CreatePDataArray(pCellsNode,
                     geosxToVTKTypeMap.at(std::type_index(typeid(integer))),
                     "types", 1, format);

    // Find all the cell fields to output
    auto pCellDataNode = pUnstructureGridNode.append_child("PCellData");
    for (auto &cellField : cellFields) {
      CreatePDataArray(pCellDataNode, std::get<1>(cellField),
                       std::get<0>(cellField), std::get<2>(cellField), format);
    }

    // Declaration of the "Piece" nodes refering to the vtu files
    for (int i = 0; i < mpiSize; i++) {
      auto curPieceNode = pUnstructureGridNode.append_child("Piece");
      string fileName = std::to_string(i) + ".vtu";
      curPieceNode.append_attribute("Source") = fileName.c_str();
    }

    // Save the files
    string pvdFileName = m_baseName + ".pvd";
    m_rootFile.save_file(pvdFileName.c_str());
    pvtuFile.save_file(pvtuFileName.c_str());
  }

  string vtuFileName =
      timeStepFolderName + "/" + std::to_string(mpiRank) + ".vtu";
  CustomVTUXMLWriter vtuWriter(vtuFileName);
  vtuWriter.WriteHeader();
  vtuWriter.OpenXMLNode("VTKFile", {{"type", "UnstructuredGrid"},
                                    {"version", "0.1"},
                                    {"byte_order", "LittleEndian"}});
  vtuWriter.OpenXMLNode("UnstructuredGrid", {});

  // Declaration of the node Piece and the basic informations of the mesh
  localIndex totalNumberOfCells =
      elemManager->getNumberOfElements<CellElementSubRegion>();
  localIndex totalNumberOfSubRegion = 0;
  elemManager->forElementRegionsComplete<CellElementRegion>(
      [&](localIndex const GEOSX_UNUSED_ARG(er), auto const *const elemRegion) {
        totalNumberOfSubRegion += elemRegion->numSubRegions();
      });
  vtuWriter.OpenXMLNode(
      "Piece", {{"NumberOfPoints", std::to_string(nodeManager->size())},
                {"NumberOfCells", std::to_string(totalNumberOfCells)}});

  // Definition of node Points
  vtuWriter.OpenXMLNode("Points", {});

  // Definition of the node DataArray that will contain all the node coordinates
  vtuWriter.OpenXMLNode(
      "DataArray",
      {{"type", geosxToVTKTypeMap.at(std::type_index(typeid(real64)))},
       {"Name", "Position"},
       {"NumberOfComponents", "3"},
       {"format", format}});
  vtuWriter.WriteVertices(nodeManager->referencePosition(), m_binary);
  vtuWriter.CloseXMLNode("DataArray");
  vtuWriter.CloseXMLNode("Points");

  // Point data output
  vtuWriter.OpenXMLNode("PointData", {});
  for (auto &nodeField : nodeFields) {
    WrapperBase const *const wrapper =
        nodeManager->getWrapperBase(std::get<0>(nodeField));
    vtuWriter.OpenXMLNode(
        "DataArray",
        {{"type", std::get<1>(nodeField)},
         {"Name", std::get<0>(nodeField)},
         {"NumberOfComponents", std::to_string(std::get<2>(nodeField))},
         {"format", format}});
    rtTypes::ApplyArrayTypeLambda1(
        std::get<3>(nodeField), [&](auto type) -> void {
          using cType = decltype(type);
          const Wrapper<cType> &view = Wrapper<cType>::cast(*wrapper);
          vtuWriter.WriteNodeData(view, m_binary);
        });
    vtuWriter.CloseXMLNode("DataArray");
  }
  vtuWriter.CloseXMLNode("PointData");

  // Definition of the node Cells
  vtuWriter.OpenXMLNode("Cells", {});

  // Definition of the node DataArray that will contain the connectivities
  vtuWriter.OpenXMLNode(
      "DataArray",
      {{"type", geosxToVTKTypeMap.at(std::type_index(typeid(localIndex)))},
       {"Name", "connectivity"},
       {"NumberOfComponents", "1"},
       {"format", format}});

  vtuWriter.WriteCellConnectivities(elemManager, m_binary);

  /*
  elemManager->forElementRegionsComplete< ElementRegion >( [&]( localIndex const
  er, auto const * const elemRegion )
  {
    elemRegion->template forElementSubRegions< CellElementSubRegion >( [&]( auto
  const * const elemSubRegion )
    {
      vtuWriter.WriteCellConnectivities( elemSubRegion->GetElementTypeString(),
  elemSubRegion->nodeList(), m_binary );
    });
  });
  */
  vtuWriter.CloseXMLNode("DataArray");

  array1d<std::tuple<integer, localIndex, string>>
      subRegionsInfo; // First value : cell size, Second value : number of
                      // cells, Third value : cell Types
  elemManager->forElementRegionsComplete<CellElementRegion>(
      [&](localIndex const GEOSX_UNUSED_ARG(er), auto const *const elemRegion) {
        elemRegion->template forElementSubRegions<CellElementSubRegion>(
            [&](auto const *const elemSubRegion) {
              subRegionsInfo.push_back(std::make_tuple(
                  elemSubRegion->numNodesPerElement(), elemSubRegion->size(),
                  elemSubRegion->GetElementTypeString()));
            });
      });

  // Definition of the node DataArray that will contain the offsets
  vtuWriter.OpenXMLNode(
      "DataArray",
      {{"type", geosxToVTKTypeMap.at(std::type_index(typeid(localIndex)))},
       {"Name", "offsets"},
       {"NumberOfComponents", "1"},
       {"format", format}});
  vtuWriter.WriteCellOffsets(elemManager, m_binary);
  vtuWriter.CloseXMLNode("DataArray");

  // Definition of the node DataArray that will contain the cell types
  vtuWriter.OpenXMLNode(
      "DataArray",
      {{"type", geosxToVTKTypeMap.at(std::type_index(typeid(integer)))},
       {"Name", "types"},
       {"NumberOfComponents", "1"},
       {"format", format}});
  vtuWriter.WriteCellTypes(elemManager, m_binary);
  vtuWriter.CloseXMLNode("DataArray");

  vtuWriter.CloseXMLNode("Cells");
  // Definition of the CellDataArray node that will contains all the data held
  // by the elements
  vtuWriter.OpenXMLNode("CellData", {});
  for (auto &cellField : cellFields) {
    vtuWriter.OpenXMLNode(
        "DataArray",
        {{"type", std::get<1>(cellField)},
         {"Name", std::get<0>(cellField)},
         {"NumberOfComponents", std::to_string(std::get<2>(cellField))},
         {"format", format}});
    rtTypes::ApplyArrayTypeLambda1(
        std::get<3>(cellField), [&](auto type) -> void {
          using cType = decltype(type);
          auto dataView =
              elemManager->ConstructViewAccessor<cType>(std::get<0>(cellField));
          vtuWriter.WriteCellData(dataView, elemManager, m_binary);
        });
    vtuWriter.CloseXMLNode("DataArray");
  }
  vtuWriter.CloseXMLNode("CellData");
  vtuWriter.CloseXMLNode("Piece");
  vtuWriter.CloseXMLNode("UnstructuredGrid");
  vtuWriter.CloseXMLNode("VTKFile");
}

} // namespace geosx
