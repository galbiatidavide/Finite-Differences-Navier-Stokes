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
    MeshHandler(Params given_input, std::string file) : GridType(given_input), fileName(file) {};

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

    double x_min = std::numeric_limits<double>::max();
    double y_min = std::numeric_limits<double>::max();
    double z_min = std::numeric_limits<double>::max();

    double x_max = std::numeric_limits<double>::lowest();
    double y_max = std::numeric_limits<double>::lowest();
    double z_max = std::numeric_limits<double>::lowest();

    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++) {
        double p[3];
        points->GetPoint(i, p);
        vertices[i] = {{ p[0], p[1], p[2] }};

        x_min = std::min(x_min, p[0]);
        y_min = std::min(y_min, p[1]);
        z_min = std::min(z_min, p[2]);

        x_max = std::max(x_max, p[0]);
        y_max = std::max(y_max, p[1]);
        z_max = std::max(z_max, p[2]);
    }

    // Calculate the maximum dimension of the bounding box
    double max_dimension = std::max({x_max - x_min, y_max - y_min, z_max - z_min});

    // Normalize vertices using the maximum dimension
    for (auto& v : vertices) {
        v[0] = (v[0] - x_min) / max_dimension;
        v[1] = (v[1] - y_min) / max_dimension;
        v[2] = (v[2] - z_min) / max_dimension;
    }

    // Extract faces
    vtkSmartPointer<vtkCellArray> triangles = polyData->GetPolys();
    triangles->InitTraversal();
    vtkIdType npts;
    const vtkIdType *pts;
    while (triangles->GetNextCell(npts, pts)) {
        if (npts == 3) {
            faces.push_back({{ pts[0], pts[1], pts[2] }});
        }
    }

    };

    void save_vertices() {
        std::cout << "Saving vertices to file" << std::endl;
        std::cout << vertices.size() << std::endl;
        std::ofstream file;
        file.open("vertices.txt");
        for (int i = 0; i < vertices.size(); i++) {
            file << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << std::endl;
        }
        file.close();
    };

    void save_faces() {
        std::ofstream facesFile("faces.txt");
        if (facesFile.is_open()) {
            for (const auto& f : faces) {
                facesFile << f[0] << " " << f[1] << " " << f[2] << std::endl;
            }
            facesFile.close();
            std::cout << "Faces saved to faces.txt" << std::endl;
        } else {
            std::cerr << "Unable to open file faces.txt for writing" << std::endl;
        }
    };


    bool isPointInsideMesh(const std::array<double, 3>& point,
                       const std::vector<std::array<double, 3>>& vertices,
                       const std::vector<std::array<int, 3>>& faces) {
    // Define a ray direction. We'll use the positive x-direction.
    std::array<double, 3> rayDirection = {1.0, 0.0, 0.0};

    int intersectionCount = 0;

    // Iterate over all faces (triangles)
    for (const auto& face : faces) {
        const std::array<double, 3>& v0 = vertices[face[0]];
        const std::array<double, 3>& v1 = vertices[face[1]];
        const std::array<double, 3>& v2 = vertices[face[2]];

        if (rayIntersectsTriangle(point, rayDirection, v0, v1, v2)) {
            intersectionCount++;
        }
    }

    };

    bool rayIntersectsTriangle(const std::array<double, 3>& rayOrigin,
                           const std::array<double, 3>& rayVector,
                           const std::array<double, 3>& v0,
                           const std::array<double, 3>& v1,
                           const std::array<double, 3>& v2) {
    const double EPSILON = 1e-8;
    std::array<double, 3> edge1, edge2, h, s, q;
    double a, f, u, v;

    // Calculate edges
    for (int i = 0; i < 3; ++i) {
        edge1[i] = v1[i] - v0[i];
        edge2[i] = v2[i] - v0[i];
    }

    // Calculate determinant
    h[0] = rayVector[1] * edge2[2] - rayVector[2] * edge2[1];
    h[1] = rayVector[2] * edge2[0] - rayVector[0] * edge2[2];
    h[2] = rayVector[0] * edge2[1] - rayVector[1] * edge2[0];

    a = edge1[0] * h[0] + edge1[1] * h[1] + edge1[2] * h[2];

    if (a > -EPSILON && a < EPSILON)
        return false; // This means the ray is parallel to the triangle.

    f = 1.0 / a;

    // Calculate u parameter and test bound
    for (int i = 0; i < 3; ++i) {
        s[i] = rayOrigin[i] - v0[i];
    }

    u = f * (s[0] * h[0] + s[1] * h[1] + s[2] * h[2]);
    if (u < 0.0 || u > 1.0)
        return false;

    // Calculate v parameter and test bound
    q[0] = s[1] * edge1[2] - s[2] * edge1[1];
    q[1] = s[2] * edge1[0] - s[0] * edge1[2];
    q[2] = s[0] * edge1[1] - s[1] * edge1[0];

    v = f * (rayVector[0] * q[0] + rayVector[1] * q[1] + rayVector[2] * q[2]);
    if (v < 0.0 || u + v > 1.0)
        return false;

    // At this stage, we can compute t to find out where the intersection point is on the line.
    double t = f * (edge2[0] * q[0] + edge2[1] * q[1] + edge2[2] * q[2]);
    if (t > EPSILON) // ray intersection
        return true;

    return false; 
};



    // Destructor
    ~MeshHandler() {};
};
