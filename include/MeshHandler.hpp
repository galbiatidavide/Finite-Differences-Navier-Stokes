#include "Grid.hpp"
#include <vtkSTLReader.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>
#include <vtkPoints.h>
#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <unordered_map>
#include <cmath>
#include <petscdm.h>
#include <petscdmstag.h>
#include <petscksp.h>



template <typename GridType>
class MeshHandler : public GridType {

protected:
std::string fileName;
std::vector<std::array<double, 3>> vertices;
std::vector<std::array<int, 3>> faces;

public:
    MeshHandler(Params given_input, std::string file) : GridType(given_input), fileName(file) {}

    // Method to read STL file
    void reader() {
        
        vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
        reader->SetFileName(fileName.c_str());
        reader->Update();

        // Get the polydata object
        vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();

        // Extract vertices
        vtkSmartPointer<vtkPoints> points = polyData->GetPoints();
        vertices.resize(points->GetNumberOfPoints());
        vtkSmartPointer<vtkCellArray> triangles = polyData->GetPolys();
        triangles->InitTraversal();

    }

    void save_vertices() {
        std::ofstream file;
        file.open("vertices.txt");
        for (int i = 0; i < vertices.size(); i++) {
            file << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << std::endl;
        }
        file.close();
    }


    // Destructor
    ~MeshHandler() {}
};
