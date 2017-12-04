/*
 * MeshManager.cpp
 *
 *  Created on: Oct 18, 2017
 *      Author: sherman
 */

#include "MeshManager.hpp"

#include "../../../cxx-utilities/src/src/DocumentationNode.hpp"
#include "MeshGeneratorBase.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

MeshManager::MeshManager( std::string const & name,
                          ManagedGroup * const parent ):
  ManagedGroup( name, parent)
{}

MeshManager::~MeshManager()
{}


void MeshManager::FillDocumentationNode( dataRepository::ManagedGroup * const /*group*/ )
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName("Mesh");
  docNode->setSchemaType("UniqueNode");
  docNode->setShortDescription("Mesh manager");

}


void MeshManager::CreateChild( string const & childKey, string const & childName )
{
  std::cout << "Adding Mesh: " << childKey << ", " << childName << std::endl;
  std::unique_ptr<MeshGeneratorBase> solver = MeshGeneratorBase::CatalogInterface::Factory( childKey, childName, this );
  this->RegisterGroup<MeshGeneratorBase>( childName, std::move(solver) );
}


void MeshManager::GenerateMeshes( dataRepository::ManagedGroup * const domain )
{
  this->forSubGroups<MeshGeneratorBase>([this, domain]( MeshGeneratorBase * meshGen ) -> void
    {
      meshGen->GenerateMesh( domain );
    });
}



} /* namespace geosx */
