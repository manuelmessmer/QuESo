
#include "geometries/element.h"

void Element::pSetSurfaceMesh(SurfaceMeshPtrType& pSurfaceMesh ){
        mpSurfaceMesh = std::move(pSurfaceMesh);
        mpInsideTest = std::make_unique<InsideTestType>(mpSurfaceMesh, mLocalLowerPoint, mLocalUpperPoint);

        mSurfaceMeshSetFlag = true;
    }

    Element::SurfaceMeshPtrType Element::pGetSurfaceMesh(){
        if( !mSurfaceMeshSetFlag ){
            throw  std::runtime_error("Element :: Surface Mesh Pointer has not been set" );
        }
        return mpSurfaceMesh;
    }

    void Element::ClearSurfaceMesh(){
        mpInsideTest.reset();
        mpSurfaceMesh = nullptr;
        //mpSurfaceMesh.reset(); // TODO!!
        mSurfaceMeshSetFlag = false;
    }

    // TODO: Use inside_test.h and remove corresponding include.
    bool Element::IsPointInTrimmedDomain(PointType& rTestPoint){
        if( !mSurfaceMeshSetFlag ){
            throw  std::runtime_error("Element :: Surface Mesh Pointer has not been set" );
        }
        std::array<double,3> tmp_point = {rTestPoint[0], rTestPoint[1], rTestPoint[2]};

        return mpInsideTest->is_inside(tmp_point);
    }

    Element::BoundingBox Element::ComputeTrimmedBoundingBox(){
        if( !mSurfaceMeshSetFlag ){
            throw  std::runtime_error("Element :: Surface Mesh Pointer has not been set" );
        }
        BoundingBox bounding_box = { {1e10, 1e10, 1e10},
                                 {-1e10, -1e10, -1e10} };

        vtkSmartPointer<vtkPoints> points = mpSurfaceMesh->GetPoints();
        vtkSmartPointer<vtkDataArray> dataArray = points->GetData();
        const int size_a = points->GetNumberOfPoints();
        const int size_b = dataArray->GetNumberOfComponents();

        for( int i = 0; i < size_a; i++){
            std::array<double, 3> tmp_coordinate= { dataArray->GetComponent(i,0),
                                                    dataArray->GetComponent(i,1),
                                                    dataArray->GetComponent(i,2)};
            for( IndexType i = 0; i < 3; ++i){
                if( tmp_coordinate[i] < bounding_box[0][i] ){ // Find min values
                    bounding_box[0][i] = tmp_coordinate[i];
                }
                if( tmp_coordinate[i] > bounding_box[1][i] ){ // Find max values
                    bounding_box[1][i] = tmp_coordinate[i];
                }
            }
        }

        return bounding_box;
    }

    Element::TriangleVectorType& Element::GetTriangles(){

        if( mTriangles.size() > 0)
            return mTriangles;

        if( !mSurfaceMeshSetFlag ){
            throw  std::runtime_error("Element :: Surface Mesh Pointer has not been set" );
        }

        // Isotropic remeshing
        auto remeshed_intersection = MeshUtilities::IsotropicRemeshing(mpSurfaceMesh, mParameters, this->GetId());
        //auto remeshed_intersection = mpSurfaceMesh;


        mTriangles.reserve(remeshed_intersection->GetNumberOfPolys());

        vtkSmartPointer<vtkPolyDataNormals> dataset = vtkSmartPointer<vtkPolyDataNormals>::New();
        dataset->SetInputData(remeshed_intersection);
        dataset->ComputePointNormalsOff();
        dataset->ComputeCellNormalsOn();
        dataset->Update();

        remeshed_intersection = dataset->GetOutput();

        auto normals = remeshed_intersection->GetCellData()->GetNormals();
        auto polys = remeshed_intersection->GetPolys();
        auto points = remeshed_intersection->GetPoints();

        const int num_cells = remeshed_intersection->GetNumberOfPolys();

        polys->InitTraversal();
        for( int i = 0; i < num_cells; i++){
            vtkIdType npts, *pts;
            auto h = polys->GetNextCell(npts, pts);
            std::array<double,3> normal{};
            normals->GetTuple(i, normal.data());

            std::array<PointType, 3> coordinates{};
            if( h == 0){
                break;
            }
            if( npts == 3 ){
                std::vector<size_t> tmp_vec = { (size_t)pts[0], (size_t)pts[1], (size_t)pts[2]};

                for( int j = 0; j < npts; j++){
                    //auto test = points->GetData()->GetTuple(pts[0]);
                    auto tuple = points->GetData()->GetTuple(pts[j]);
                    coordinates[j] = {tuple[0], tuple[1], tuple[2]};
                }

                // Construct triangle to get integration points
                Triangle3D3N tmp_triangle(coordinates[0], coordinates[1], coordinates[2], normal);
                if( tmp_triangle.Area() > 1e-12 )
                    mTriangles.push_back(tmp_triangle);
            }
        }
        //this->ClearSurfaceMesh();
        return mTriangles;
    }
