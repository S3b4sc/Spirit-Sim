#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <map>
#include <array>
#include <cmath>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

using Kernel       = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point2       = Kernel::Point_2;
using Triangulation = CGAL::Delaunay_triangulation_2<Kernel>;



// ----- Data Structures -----
struct Vec3 {
    double x, y, z;
    Vec3 operator-(const Vec3& o) const { 
        return {x-o.x, y-o.y, z-o.z}; 
    }

    Vec3 cross(const Vec3& o) const {
        return { y*o.z - z*o.y,
                 z*o.x - x*o.z,
                 x*o.y - y*o.x };
    }

    double dot(const Vec3& o) const { 
        return x*o.x + y*o.y + z*o.z; 
    }

    Vec3 normalize() const {
        double n = std::sqrt(x*x + y*y + z*z);
        return {x/n, y/n, z/n};
    }

};


// Triangle struct for contributions
struct Triangle { int i0, i1, i2; };

// Generate Spirit-style triangles for a 2D lattice
// positions: original atom positions (size = n_basis_atoms)
// n_basis_atoms: number of basis atoms per cell
// Lx, Ly: number of cells in x and y
// pbc_x, pbc_y: periodic boundary conditions flags
std::vector<Triangle> generate_spirit_triangles(
    const std::vector<Vec3>& positions,
    int n_basis_atoms,
    int Lx, int Ly,
    bool pbc_x, bool pbc_y)
{
    // Stretch factor for boundary anchors
    double stretch = 0.1;
    // compute lattice vectors (assuming unit vectors along axes)
    Vec3 ta{1,0,0}, tb{0,1,0};

    // Build basis cell points: basis atoms + 3 corners
    int Nbase = n_basis_atoms + 3;
    std::vector<Vec3> base_pts(Nbase);
    for(int i=0;i<n_basis_atoms;++i) base_pts[i] = {positions[i].x, positions[i].y, 0.0};
    // corner a+b
    Vec3 bb = { ta.x+tb.x, ta.y+tb.y, 0.0 };
    base_pts[n_basis_atoms]   = {positions[0].x+bb.x*stretch, positions[0].y+bb.y*stretch, 0.0};
    // corner b
    base_pts[n_basis_atoms+1] = {positions[0].x+tb.x*(1-stretch), positions[0].y+tb.y*(1-stretch),0.0};
    // corner a
    base_pts[n_basis_atoms+2] = {positions[0].x+ta.x*(1-stretch), positions[0].y+ta.y*(1-stretch),0.0};

    // Compute Delaunay triangulation on base_pts with CGAL
    using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Vb = CGAL::Triangulation_vertex_base_with_info_2<std::size_t, Kernel>;
    using Fb = CGAL::Triangulation_face_base_2<Kernel>;
    using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
    using Triangulation = CGAL::Delaunay_triangulation_2<Kernel, Tds>;
    using Point2 = Kernel::Point_2;
    
    // prepare CGAL points with info
    std::vector<std::pair<Point2,std::size_t>> pts;
    pts.reserve(Nbase);
    for(size_t i=0;i<base_pts.size();++i)
        pts.emplace_back(Point2(base_pts[i].x, base_pts[i].y), i);
    Triangulation dt; dt.insert(pts.begin(), pts.end());

    // collect base triangles (info holds base index)
    std::vector<std::array<int,3>> base_tri;
    for(auto fit=dt.finite_faces_begin(); fit!=dt.finite_faces_end(); ++fit){
        base_tri.push_back({{
            int(fit->vertex(0)->info()),
            int(fit->vertex(1)->info()),
            int(fit->vertex(2)->info())
        }});
    }

    // Now translate base triangles into global lattice
    std::vector<Triangle> all_tri;
    all_tri.reserve(base_tri.size() * Lx * Ly);
    for(int b=0;b<Ly;++b){
        for(int a=0;a<Lx;++a){
            bool a_ok = (a+1< Lx) || pbc_x;
            bool b_ok = (b+1< Ly) || pbc_y;
            for(auto &tri : base_tri){
                bool valid=true;
                std::array<int,3> idx;
                for(int i=0;i<3;++i){
                    int vid = tri[i];
                    if(vid < n_basis_atoms){
                        idx[i] = vid + a*n_basis_atoms + b*n_basis_atoms*Lx;
                    } else if(vid==n_basis_atoms && a_ok && b_ok){
                        // a+b
                        idx[i] = ((a+1)%Lx)*n_basis_atoms + ((b+1)%Ly)*n_basis_atoms*Lx;
                    } else if(vid==n_basis_atoms+1 && b_ok){
                        // b
                        idx[i] = a*n_basis_atoms + ((b+1)%Ly)*n_basis_atoms*Lx;
                    } else if(vid==n_basis_atoms+2 && a_ok){
                        // a
                        idx[i] = ((a+1)%Lx)*n_basis_atoms + b*n_basis_atoms*Lx;
                    } else {
                        valid=false; break;
                    }
                }
                if(valid) all_tri.push_back({idx[0],idx[1],idx[2]});
            }
        }
    }
    return all_tri;
}


// --- Solid angle (Berg–Lüscher) ---
double solid_angle(const Vec3& S1, const Vec3& S2, const Vec3& S3) {
    Vec3 a = S1.normalize();
    Vec3 b = S2.normalize();
    Vec3 c = S3.normalize();

    // Compute the solid angle using the Berg–Lüscher formula
    double triple = a.dot(b.cross(c));
    double denom  = 1.0 + a.dot(b) + b.dot(c) + c.dot(a);
    return 2.0 * std::atan2(triple, denom);
}


// ----- Load Data -----
bool load_spins(const std::string& filename, std::vector<Vec3>& spins) {
    
    std::ifstream infile(filename);

    if (!infile) return false;
    std::string line;

    while (std::getline(infile, line)) {
        // Skip empty lines or comment lines
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        double sx, sy, sz;
        if (ss >> sx >> sy >> sz)
            spins.push_back({sx, sy, sz});
    }
    return true;
}

bool load_positions(const std::string& filename,
                std::vector<Vec3>& positions)
{
    std::ifstream infile(filename);
    if (!infile) return false;

    std::string line;
    while (std::getline(infile, line)) {

        // Skip empty lines or comment lines
        if (line.empty() || line[0] == '#') continue;

        // Parse the line into position components 
        std::stringstream ss(line);
        double x, y, z;
        if (ss >> x >> y >> z) {
            positions.push_back({x, y, z});
        }
    }

    return true;
}


// Slice layers by z value (within tolerance)
std::map<double, std::vector<std::pair<Vec3, Vec3>>> split_by_layers(
    const std::vector<Vec3>& positions, const std::vector<Vec3>& spins,
    double z_tol = 1e-3)
{
    std::map<double, std::vector<std::pair<Vec3, Vec3>>> layers;

    for (size_t i = 0; i < positions.size(); ++i) {
        double z = positions[i].z;

        // Round z to nearest tolerance
        double z_key = std::round(z / z_tol) * z_tol;
        layers[z_key].emplace_back(positions[i], spins[i]);
    }

    return layers;
}

double compute_topological_charge_2D(
    const std::vector<Vec3>& positions2D,
    const std::vector<Vec3>& spins2D)
{
    //Triangulation dt;
    //std::vector<Point2> points2D;
//
    //// Convert Vec3 positions to Point2 for 2D triangulation
    //for (const auto& p : positions2D)
    //    points2D.emplace_back(p.x, p.y);
//
    //// Insert points into the Delaunay triangulation
    //dt.insert(points2D.begin(), points2D.end());

    //double Q = 0.0;
    //const double four_pi = 4.0 * M_PI;
//
    //for (auto fit = dt.finite_faces_begin(); fit != dt.finite_faces_end(); ++fit) {
    //    Point2 p0 = fit->vertex(0)->point();
    //    Point2 p1 = fit->vertex(1)->point();
    //    Point2 p2 = fit->vertex(2)->point();
//
    //    auto find_index = [&](const Point2& p) {
    //        // Find the index of the point in positions2D
    //        // using a simple linear search
    //        for (size_t i = 0; i < positions2D.size(); ++i)
    //            if (std::abs(positions2D[i].x - p.x()) < 1e-6 &&
    //                std::abs(positions2D[i].y - p.y()) < 1e-6)
    //                return int(i);
    //        return -1;
    //    };
//
    //    int i0 = find_index(p0), i1 = find_index(p1), i2 = find_index(p2);
    //    if (i0 < 0 || i1 < 0 || i2 < 0) continue;
//
        // Compute the signed area of the triangle
        // using the cross product of vectors on 2D
        // (p1 - p0) and (p2 - p0)
        // This gives us the orientation and area, CW is negative, CCW is positive
        //double zx = (p1.x() - p0.x()) * (p2.y() - p0.y()) -
        //            (p1.y() - p0.y()) * (p2.x() - p0.x());
        //double zx = ((p0.x() - p1.x()) * (p0.y() - p2.y())) - ((p0.y() - p1.y()) * (p0.x() - p2.x()));
        //double sign = (zx > 0 ? +1.0 : -1.0);
        //Vec3 v0 = { p0.x(), p0.y(), 0.0 };
        //Vec3 v1 = { p1.x(), p1.y(), 0.0 };
        //Vec3 v2 = { p2.x(), p2.y(), 0.0 };
//
        //Vec3 n = (v0 - v1).cross(v0 - v2).normalize();
        //// The sign of the area is determined by the z-component of the normal vector
        //double sign = (std::abs(n.z) > 1e-10) ? (n.z / std::abs(n.z)) : 0.0;
//
        //double omega = solid_angle(spins2D[i0], spins2D[i1], spins2D[i2]);
        //std::cout << "Triangle: " << i0 << ", " << i1 << ", " << i2
        //          << " → Area sign: " << sign
        //          << ", Solid angle: " << omega << "\n";
        //Q += sign * (omega / four_pi);
    //}
    std::vector<Triangle> triangles = generate_spirit_triangles(
    positions,      // base cell atom positions (size = n_basis_atoms)
    n_basis_atoms,
    Lx, Ly,
    pbc_x, pbc_y
    );
    double Q = 0.0;
    for (const auto& tri : triangles) {
        int i0 = tri.i0;
        int i1 = tri.i1;
        int i2 = tri.i2;
    
        // Use spins[i0], spins[i1], spins[i2]
        double omega = solid_angle(spins[i0], spins[i1], spins[i2]);
    
        // Use Spirit's sign convention if needed
        double dQ = omega / (4.0 * M_PI);
        Q += dQ;
    }

    return Q;
}

// Entry point to call from Python
extern "C"
void compute_layered_topo_charge(const char* pos_file, const char* spin_file) {
    std::vector<Vec3> positions, spins;

    if (!load_positions(pos_file, positions)) {
        std::cerr << "Failed to load positions\n";
        return;
    }
    if (!load_spins(spin_file, spins)) {
        std::cerr << "Failed to load spins\n";
        return;
    }

    if (positions.size() != spins.size()) {
        std::cerr << "Mismatched number of positions and spins\n";
        return;
    }

    auto layers = split_by_layers(positions, spins);

    for (const auto& [z, layer] : layers) {

        std::vector<Vec3> pos2D, spin2D;

        for (const auto& [p, s] : layer) {
            pos2D.push_back(p);
            spin2D.push_back(s);
        }

        double Q = compute_topological_charge_2D(pos2D, spin2D);
        std::cout << "Layer z=" << z << " → Q = " << Q << "\n";
    }
}

// Compile using: g++ -O3 -std=c++17 -fPIC -shared topological_charge.cpp -o libtopo.so -lgmp -lmpfr


