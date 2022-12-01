// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de


//// Project includes
#include "io/io_utilities.h"

namespace tibra {

bool IO::WriteMeshToSTL(const TriangleMesh& rTriangleMesh,
                        const char* Filename,
                        const bool Binary){
  std::ofstream file;
  if(Binary)
    file.open(Filename, std::ios::out | std::ios::binary);
  else
    file.open(Filename);

  if(!file.good()){
    std::cerr << "Warning :: IO::WriteMeshToSTL :: Could not open file: " << Filename << ".\n";
    return false;
  }

  const IndexType num_triangles = rTriangleMesh.NumOfTriangles();
  if(Binary)
  {
    file << "FileType: Binary                                                                ";
    file.write(reinterpret_cast<const char *>(&num_triangles), sizeof(num_triangles));

    for(int triangle_id = 0; triangle_id < num_triangles; ++triangle_id)
    {
      const auto& p1 = rTriangleMesh.P1(triangle_id);
      const auto& p2 = rTriangleMesh.P2(triangle_id);
      const auto& p3 = rTriangleMesh.P3(triangle_id);

      const auto& normal = rTriangleMesh.Normal(triangle_id);

      const float coords[12] = { static_cast<float>(normal[0]), static_cast<float>(normal[1]), static_cast<float>(normal[2]),
                                 static_cast<float>(p1[0]), static_cast<float>(p1[1]), static_cast<float>(p1[2]),
                                 static_cast<float>(p2[0]), static_cast<float>(p2[1]), static_cast<float>(p2[2]),
                                 static_cast<float>(p3[0]), static_cast<float>(p3[1]), static_cast<float>(p3[2]) };

      for(int i=0; i<12; ++i)
        file.write(reinterpret_cast<const char *>(&coords[i]), sizeof(coords[i]));
      file << "  ";
    }
  }
  else
  {
    file << "solid\n";
    for(int triangle_id = 0; triangle_id < num_triangles; ++triangle_id)
    {
      const auto& p1 = rTriangleMesh.P1(triangle_id);
      const auto& p2 = rTriangleMesh.P2(triangle_id);
      const auto& p3 = rTriangleMesh.P3(triangle_id);

      const auto& normal = rTriangleMesh.Normal(triangle_id);

      file << "facet normal " << normal[0] << ' ' << normal[1] << ' ' << normal[2] << "\nouter loop\n";
      file << "vertex " << p1[0] << ' ' << p1[1] << ' ' << p1[2] << ' ' << "\n";
      file << "vertex " << p2[0] << ' ' << p2[1] << ' ' << p2[2] << ' ' << "\n";
      file << "vertex " << p3[0] << ' ' << p3[1] << ' ' << p3[2] << ' ' << "\n";
      file << "endloop\nendfacet\n";
    }
    file << "endsolid"<<std::endl;
  }

  return true;
}

bool IO::ReadMeshFromSTL(TriangleMesh& rTriangleMesh,
                         const char* Filename){

  // Open file
  std::ifstream file(Filename, std::ios::binary);

  if( !file.good() ) {
    std::cerr << "Warning :: IO::ReadMeshFromSTL :: Couldnt handle file: " << Filename << ". Please provide .stl as binary.\n";
    return false;
  }

  // Ignore the first 80 chars of the header
  int position = 0;
  char message;
  std::string test_binary_ascii{};
  while(position < 80) {
    file.read(reinterpret_cast<char*>(&message), sizeof(message));
    if( position < 5 ){
      test_binary_ascii.push_back(message);
    }

    if(!file.good())
      break;

    ++position;
  }

  if( test_binary_ascii == "solid" ) { // If the first 5 characters are "solid"
    std::cerr << "Warning :: IO::ReadMeshFromSTL :: Read STL from ASCII is not implemented yet. Please use binary format.\n";
    return false;
  }

  if(position != 80) {
    std::cerr << "Warning :: IO::ReadMeshFromSTL :: File " << Filename << " is empty.\n";
    return false;
  }

  int index = 0;
  std::map<Vector3d, IndexType> index_map;

  // Read number of triangles
  unsigned int num_triangles;
  if(!(file.read(reinterpret_cast<char*>(&num_triangles), sizeof(num_triangles)))) {
    std::cerr << "Warning :: IO::ReadMeshFromSTL :: Couldnt read number of triangles. \n";
    return false;
  }
  rTriangleMesh.Clear();
  rTriangleMesh.Reserve(num_triangles);

  // Loop over all triangles
  for(unsigned int i=0; i<num_triangles; ++i) {
    // Read normals
    float normal[3];
    if(!(file.read(reinterpret_cast<char*>(&normal[0]), sizeof(normal[0]))) ||
        !(file.read(reinterpret_cast<char*>(&normal[1]), sizeof(normal[1]))) ||
        !(file.read(reinterpret_cast<char*>(&normal[2]), sizeof(normal[2])))) {
      std::cerr << "Warning :: IO::ReadMeshFromSTL :: Couldnt read normals. \n";
      return false;
    }
    rTriangleMesh.AddNormal( {normal[0], normal[1], normal[2]} );

    // Read triangles and vertices. Each vertex is read seperately.
    Vector3i triangle{};
    for(int j=0; j<3; ++j) {
      float x,y,z;
      if(!(file.read(reinterpret_cast<char*>(&x), sizeof(x))) ||
          !(file.read(reinterpret_cast<char*>(&y), sizeof(y))) ||
          !(file.read(reinterpret_cast<char*>(&z), sizeof(z)))) {
        std::cerr << "Warning :: IO::ReadMeshFromSTL :: Couldnt read coordinates. \n";
        return false;
      }

      Vector3d vertex = {x,y,z};

      // Map is used to ensure unique vertices. Note that STL does not reuse vertices.
      auto index_map_iterator = index_map.insert(std::make_pair(vertex, -1)).first;
      if(index_map_iterator->second == -1) {
        triangle[j] = index;
        index_map_iterator->second = index++;
        rTriangleMesh.AddVertex(vertex);
      }
      else {
        triangle[j] = index_map_iterator->second;
      }
    }
    rTriangleMesh.AddTriangle(triangle);

    // Read so-called attribute byte count and ignore it
    char c;
    if(!(file.read(reinterpret_cast<char*>(&c), sizeof(c))) ||
        !(file.read(reinterpret_cast<char*>(&c), sizeof(c)))) {
      std::cerr << "Warning :: IO::ReadMeshFromSTL :: Couldnt read attribute byte count.\n";
      return false;
    }
  }
  return rTriangleMesh.Check();
}

bool IO::WriteMeshToVTK(const TriangleMesh& rTriangleMesh,
                        const char* Filename,
                        const bool Binary)
{

  const SizeType num_elements = rTriangleMesh.NumOfTriangles();
  const SizeType num_points = rTriangleMesh.NumOfVertices();

  std::ofstream file;
  if(Binary)
    file.open(Filename, std::ios::out | std::ios::binary);
  else
    file.open(Filename);

  file << "# vtk DataFile Version 4.1" << std::endl;
  file << "vtk output" << std::endl;
  if(Binary)
    file << "BINARY"<< std::endl;
  else
    file << "ASCII"<< std::endl;

  file << "DATASET UNSTRUCTURED_GRID" << std::endl;
  file << "POINTS " << num_points << " double" << std::endl;


  IndexType inum = 0;
  const auto& r_vertices = rTriangleMesh.GetVertices();
  const auto v_it_begin = r_vertices.begin();
  for(IndexType i = 0; i < num_points; ++i)
  {
    const auto v_it = v_it_begin+i;
    if( Binary ){
      double rx = (*v_it)[0];
      double ry = (*v_it)[1];
      double rz = (*v_it)[2];

      WriteBinary(file, rx);
      WriteBinary(file, ry);
      WriteBinary(file, rz);
    }
    else {
      file << (*v_it)[0] << ' ' << (*v_it)[1] << ' ' << (*v_it)[2] << std::endl;
    }
  }
  file << std::endl;

  // Write Cells
  file << "Cells " << num_elements << " " << num_elements*4 << std::endl;

  for(IndexType i = 0; i < num_elements; ++i)
  {
    if( Binary ){
      int k = 3;
      WriteBinary(file, k);

      for( auto id : rTriangleMesh.VertexIds(i) ){
        k = id;
        WriteBinary(file, k);
      }
    }
    else {
      file << 3;
      for( auto id : rTriangleMesh.VertexIds(i) ){
          file << ' ' << id;
      }
      file << std::endl;
    }
  }
  file << std::endl;

  file << "CELL_TYPES " << num_elements << std::endl;
  for( int i = 0; i < num_elements; ++i){
    if( Binary ){
        int k = 5;
        WriteBinary(file, k);
    }
    else {
      file << 5 << std::endl;
    }
  }
  file << std::endl;
  file.close();

  return true;
}

bool IO::WriteDisplacementToVTK(const std::vector<Vector3d>& rDisplacement,
                                const char* Filename,
                                const bool Binary){

  const SizeType num_points = rDisplacement.size();

  std::ofstream file;
  if(Binary)
    file.open(Filename, std::ios::app | std::ios::binary);
  else
    file.open(Filename);

  file << "POINT_DATA " << num_points << std::endl;
  file << "VECTORS Displacement double" << std::endl;
  for(int i = 0; i < num_points; ++i){
      if( Binary ){
        double rw1 = rDisplacement[i][0];
        WriteBinary(file, rw1);
        double rw2 = rDisplacement[i][1];
        WriteBinary(file, rw2);
        double rw3 = rDisplacement[i][2];
        WriteBinary(file, rw3);
      }
      else {
        std::cerr << "Warning :: IO::DisplacementToVTK :: Ascii export not implemented yet. \n";
        return false;
      }
  }
  file << std::endl;
  file.close();

  return true;
}

bool IO::WriteElementsToVTK(const ElementContainer& rElementContainer, //PolygonMesh
                            const char* Filename,
                            const bool Binary){


  const SizeType num_elements = rElementContainer.size();

  std::ofstream file;
  if(Binary)
    file.open(Filename, std::ios::out | std::ios::binary);
  else
    file.open(Filename);

  file << "# vtk DataFile Version 4.1" << std::endl;
  file << "vtk output" << std::endl;
  if(Binary)
    file << "BINARY"<< std::endl;
  else
    file << "ASCII"<< std::endl;

  file << "DATASET UNSTRUCTURED_GRID" << std::endl;
  file << "POINTS " << num_elements*8 << " double" << std::endl;

  const auto begin_el_itr = rElementContainer.begin();
  for( int i = 0; i < rElementContainer.size(); ++i){
    auto el_itr = *(begin_el_itr + i);
    const auto& lower_point = el_itr->GetLowerBound();
    const auto& upper_point = el_itr->GetUpperBound();

    if( Binary ){

      double rx0 = lower_point[0];
      SwapEnd(rx0);
      double rx1 = upper_point[0];
      SwapEnd(rx1);
      double ry0 = lower_point[1];
      SwapEnd(ry0);
      double ry1 = upper_point[1];
      SwapEnd(ry1);
      double rz0 = lower_point[2];
      SwapEnd(rz0);
      double rz1 = upper_point[2];
      SwapEnd(rz1);

      file.write(reinterpret_cast<char*>(&rx0), sizeof(double));
      file.write(reinterpret_cast<char*>(&ry0), sizeof(double));
      file.write(reinterpret_cast<char*>(&rz0), sizeof(double));

      file.write(reinterpret_cast<char*>(&rx1), sizeof(double));
      file.write(reinterpret_cast<char*>(&ry0), sizeof(double));
      file.write(reinterpret_cast<char*>(&rz0), sizeof(double));

      file.write(reinterpret_cast<char*>(&rx1), sizeof(double));
      file.write(reinterpret_cast<char*>(&ry1), sizeof(double));
      file.write(reinterpret_cast<char*>(&rz0), sizeof(double));

      file.write(reinterpret_cast<char*>(&rx0), sizeof(double));
      file.write(reinterpret_cast<char*>(&ry1), sizeof(double));
      file.write(reinterpret_cast<char*>(&rz0), sizeof(double));

      file.write(reinterpret_cast<char*>(&rx0), sizeof(double));
      file.write(reinterpret_cast<char*>(&ry0), sizeof(double));
      file.write(reinterpret_cast<char*>(&rz1), sizeof(double));

      file.write(reinterpret_cast<char*>(&rx1), sizeof(double));
      file.write(reinterpret_cast<char*>(&ry0), sizeof(double));
      file.write(reinterpret_cast<char*>(&rz1), sizeof(double));

      file.write(reinterpret_cast<char*>(&rx1), sizeof(double));
      file.write(reinterpret_cast<char*>(&ry1), sizeof(double));
      file.write(reinterpret_cast<char*>(&rz1), sizeof(double));

      file.write(reinterpret_cast<char*>(&rx0), sizeof(double));
      file.write(reinterpret_cast<char*>(&ry1), sizeof(double));
      file.write(reinterpret_cast<char*>(&rz1), sizeof(double));
    }
    else {
      file << lower_point[0] << ' ' << lower_point[1] << ' ' << lower_point[2] << std::endl;
      file << upper_point[0] << ' ' << lower_point[1] << ' ' << lower_point[2] << std::endl;
      file << upper_point[0] << ' ' << upper_point[1] << ' ' << lower_point[2] << std::endl;
      file << lower_point[0] << ' ' << upper_point[1] << ' ' << lower_point[2] << std::endl;
      file << lower_point[0] << ' ' << lower_point[1] << ' ' << upper_point[2] << std::endl;
      file << upper_point[0] << ' ' << lower_point[1] << ' ' << upper_point[2] << std::endl;
      file << upper_point[0] << ' ' << upper_point[1] << ' ' << upper_point[2] << std::endl;
      file << lower_point[0] << ' ' << upper_point[1] << ' ' << upper_point[2] << std::endl;
    }
  }
  file << std::endl;
  // Write Cells
  file << "Cells " << num_elements << " " << num_elements*9 << std::endl;
  for( int i = 0; i < static_cast<int>(rElementContainer.size()); ++i){
    if( Binary ){
      int k = 8;
      WriteBinary(file, k);
      for( int j = 0; j < 8; ++j){
        k = 8*i+j;
        WriteBinary(file, k);
      }
    }
    else {
      file << 8 << ' ' << 8*i     << ' ' << 8*i + 1 << ' ' << 8*i + 2 << ' ' << 8*i + 3
                << ' ' << 8*i + 4 << ' ' << 8*i + 5 << ' ' << 8*i + 6 << ' ' << 8*i + 7 << std::endl;
    }
  }
  file << std::endl;

  file << "CELL_TYPES " << rElementContainer.size() << std::endl;
  for( int i = 0; i < static_cast<int>(rElementContainer.size()); ++i){
    if( Binary ){
        int k = 12;
        WriteBinary(file, k);
    }
    else {
      file << 12 << std::endl;
    }
  }
  file << std::endl;
  file.close();

  return true;
}


bool IO::WritePointsToVTK(const ElementContainer& rElementContainer,
                          const char* type,
                          const char* Filename,
                          const bool Binary){

  auto p_points = rElementContainer.pGetPoints(type);
  const auto begin_points_it_ptr = p_points->begin();
  const int num_points = p_points->size();
  const int num_elements = p_points->size();

  std::ofstream file;
  if(Binary)
    file.open(Filename, std::ios::out | std::ios::binary);
  else
    file.open(Filename);

  file << "# vtk DataFile Version 4.1" << std::endl;
  file << "vtk output" << std::endl;
  if(Binary)
    file << "BINARY"<< std::endl;
  else
    file << "ASCII"<< std::endl;


  file << "DATASET UNSTRUCTURED_GRID" << std::endl;
  file << "POINTS " << num_points << " double" << std::endl;

  const Parameters& param = (*rElementContainer.begin())->GetParameters();
  for(int i = 0; i < num_points; ++i){
    auto points_it = (begin_points_it_ptr + i);
    auto point_global = MappingUtilities::FromLocalToGlobalSpace(*points_it, param.PointA(), param.PointB() );

    if( Binary ){
      WriteBinary(file, point_global[0]);
      WriteBinary(file, point_global[1]);
      WriteBinary(file, point_global[2]);
    }
    else {
      file << point_global[0] << ' ' << point_global[1] << ' ' << point_global[2] << std::endl;
    }
  }
  file << std::endl;

  //Write Cells
  file << "Cells " << num_elements << " " << num_elements*2 << std::endl;
  for( int i = 0; i < num_elements; ++i){
    if( Binary ){
      int k = 1;
      WriteBinary(file, k);
      k = i;
      WriteBinary(file, k);
    }
    else {
      file << 1 << ' ' << i << std::endl;
    }
  }
  file << std::endl;

  file << "CELL_TYPES " << num_elements << std::endl;
  for( int i = 0; i < num_elements; ++i){
    if( Binary ){
        int k = 1;
        WriteBinary(file, k);
    }
    else {
      file << 1 << std::endl;
    }
  }
  file << std::endl;

  file << "POINT_DATA " << num_points << std::endl;
  file << "SCALARS Weights double 1" << std::endl;
  file << "LOOKUP_TABLE default" << std::endl;
  for(int i = 0; i < num_points; ++i){
      auto points_it = (begin_points_it_ptr + i);

      if( Binary ){
        double rw = points_it->GetWeight();
        WriteBinary(file, rw);
      }
      else {
        file << points_it->GetWeight() << std::endl;
      }
  }
  file << std::endl;
  file.close();

  return true;
}


bool IO::WritePointsToVTK(const std::vector<BoundaryIntegrationPoint>& rPoints,
                          const char* Filename,
                          const bool Binary){

  const auto begin_points_it_ptr = rPoints.begin();
  const int num_points = rPoints.size();
  const int num_elements = rPoints.size();

  std::ofstream file;
  if(Binary)
    file.open(Filename, std::ios::out | std::ios::binary);
  else
    file.open(Filename);

  file << "# vtk DataFile Version 4.1" << std::endl;
  file << "vtk output" << std::endl;
  if(Binary)
    file << "BINARY"<< std::endl;
  else
    file << "ASCII"<< std::endl;


  file << "DATASET UNSTRUCTURED_GRID" << std::endl;
  file << "POINTS " << num_points << " double" << std::endl;

  for(int i = 0; i < num_points; ++i){
    auto points_it = (begin_points_it_ptr + i);

    if( Binary ){
      // Make sure to create copy of points. "WriteBinary" will change them.
      auto p_x = (*points_it)[0];
      auto p_y = (*points_it)[1];
      auto p_z = (*points_it)[2];
      WriteBinary(file, p_x);
      WriteBinary(file, p_y);
      WriteBinary(file, p_z);
    }
    else {
      file << (*points_it)[0] << ' ' << (*points_it)[1] << ' ' << (*points_it)[2] << std::endl;
    }
  }
  file << std::endl;

  //Write Cells
  file << "Cells " << num_elements << " " << num_elements*2 << std::endl;
  for( int i = 0; i < num_elements; ++i){
    if( Binary ){
      int k = 1;
      WriteBinary(file, k);
      k = i;
      WriteBinary(file, k);
    }
    else {
      file << 1 << ' ' << i << std::endl;
    }
  }
  file << std::endl;

  file << "CELL_TYPES " << num_elements << std::endl;
  for( int i = 0; i < num_elements; ++i){
    if( Binary ){
        int k = 1;
        WriteBinary(file, k);
    }
    else {
      file << 1 << std::endl;
    }
  }
  file << std::endl;

  file << "POINT_DATA " << num_points << std::endl;
  file << "SCALARS Weights double 1" << std::endl;
  file << "LOOKUP_TABLE default" << std::endl;
  for(int i = 0; i < num_points; ++i){
      auto points_it = (begin_points_it_ptr + i);

      if( Binary ){
        double rw = points_it->GetWeight();
        WriteBinary(file, rw);
      }
      else {
        file << points_it->GetWeight() << std::endl;
      }
  }
  file << std::endl;
  file.close();

  return true;
}

} // End namespace tibra


