#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <map>
#include <cmath>
#include <numeric>
#include <string>
#include <algorithm>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point2 = Kernel::Point_2;
using Vb = CGAL::Triangulation_vertex_base_with_info_2<std::size_t, Kernel>;
using Fb = CGAL::Triangulation_face_base_2<Kernel>;
using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using Triangulation = CGAL::Delaunay_triangulation_2<Kernel, Tds>;

struct Vec3 { double x, y, z; };
using vector2_t = std::array<double, 2>;
using triangle_t = std::array<int, 3>;

// GPLM Solid angle
namespace Vectormath {
    struct Vector3 { 
        double x,y,z;
        Vector3 operator-(const Vector3&o)const{return{x-o.x,y-o.y,z-o.z};}
        Vector3 cross(const Vector3&o)const{return{y*o.z-z*o.y,z*o.x-x*o.z,x*o.y-y*o.x};}
        double dot(const Vector3&o)const{return x*o.x+y*o.y+z*o.z;}
        Vector3 normalize()const{double n=std::sqrt(x*x+y*y+z*z);return{ x/n,y/n,z/n };}
    };
    static double solid_angle_2(const Vector3&S1,const Vector3&S2,const Vector3&S3){
        auto a=S1.normalize(), b=S2.normalize(), c=S3.normalize();
        double triple=a.dot(b.cross(c));
        double denom=1+a.dot(b)+b.dot(c)+c.dot(a);
        return 2.0*std::atan2(triple, denom);
    }
}

// Compute Delaunay of a set of 2D points, returning base-cell triangles
static std::vector<triangle_t> compute_delaunay_triangulation_2D(const std::vector<vector2_t>& pts){
    std::vector<std::pair<Point2,std::size_t>> v;
    v.reserve(pts.size());
    for(size_t i=0;i<pts.size();++i) 
        v.emplace_back(Point2(pts[i][0],pts[i][1]), i);
    Triangulation dt; dt.insert(v.begin(),v.end());
    std::vector<triangle_t> tris;
    for(auto f=dt.finite_faces_begin();f!=dt.finite_faces_end();++f)
        tris.push_back({{
            int(f->vertex(0)->info()),
            int(f->vertex(1)->info()),
            int(f->vertex(2)->info())
        }});
    return tris;
}

// Spirit-style base-cell triangulation + tiling
static std::vector<triangle_t> generate_spirit_triangles(
    const std::vector<vector2_t>& base_pos,
    int n_basis,
    int Lx, int Ly,
    bool pbc_x, bool pbc_y,
    const vector2_t& a1, const vector2_t& a2)
{
    std::vector<vector2_t> B = base_pos;
    // Add periodic corner points
    B.push_back({ base_pos[0][0] + a1[0] + a2[0], base_pos[0][1] + a1[1] + a2[1] }); // (1,1)
    B.push_back({ base_pos[0][0] + a2[0], base_pos[0][1] + a2[1] });                  // (0,1)
    B.push_back({ base_pos[0][0] + a1[0], base_pos[0][1] + a1[1] });                  // (1,0)

    auto base_tri = compute_delaunay_triangulation_2D(B);
    std::vector<triangle_t> all;

    for(int b=0; b<Ly; ++b) {
        for(int a=0; a<Lx; ++a) {
            bool ok_a = pbc_x || (a+1 < Lx);
            bool ok_b = pbc_y || (b+1 < Ly);
            
            for(auto& t : base_tri) {
                std::array<int,3> idx;
                bool valid = true;
                
                for(int i=0; i<3; ++i) {
                    int v_idx = t[i];
                    if(v_idx < n_basis) {
                        idx[i] = v_idx + a*n_basis + b*n_basis*Lx;
                    }
                    else if(v_idx == n_basis) { // (1,1) corner
                        if(ok_a && ok_b) 
                            idx[i] = ((a+1)%Lx)*n_basis + ((b+1)%Ly)*n_basis*Lx;
                        else 
                            { valid = false; break; }
                    }
                    else if(v_idx == n_basis+1) { // (0,1) corner
                        if(ok_b) 
                            idx[i] = a*n_basis + ((b+1)%Ly)*n_basis*Lx;
                        else 
                            { valid = false; break; }
                    }
                    else if(v_idx == n_basis+2) { // (1,0) corner
                        if(ok_a) 
                            idx[i] = ((a+1)%Lx)*n_basis + b*n_basis*Lx;
                        else 
                            { valid = false; break; }
                    }
                    else {
                        valid = false;
                        break;
                    }
                }
                
                if(valid) 
                    all.push_back({idx[0], idx[1], idx[2]});
            }
        }
    }
    return all;
}

static double TopologicalCharge_FiniteLayer(
    const std::vector<Vectormath::Vector3>& spins,
    const std::vector<vector2_t>& positions,
    const std::vector<triangle_t>& tris)
{
    double Q = 0.0;
    const double inv4pi = 1.0/(4.0*M_PI);
    int skipped = 0;

    for(auto& t : tris) {
        int i0 = t[0], i1 = t[1], i2 = t[2];
        const auto& p0 = positions[i0], p1 = positions[i1], p2 = positions[i2];
        
        // Calculate triangle area in 2D for orientation
        double zx = (p1[0]-p0[0])*(p2[1]-p0[1]) - (p1[1]-p0[1])*(p2[0]-p0[0]);
        
        // Skip degenerate triangles
        if(std::abs(zx) < 1e-10) {
            skipped++;
            continue;
        }

        double sign = zx > 0 ? 1.0 : -1.0;
        double omega = Vectormath::solid_angle_2(spins[i0], spins[i1], spins[i2]);
        Q += sign * omega * inv4pi;
    }

    if(skipped > 0)
        std::cout << "Skipped " << skipped << " degenerate triangles\n";
    return Q;
}

static bool load_txt(const char* fn, std::vector<std::array<double,3>>& out){
    std::ifstream in(fn); 
    if(!in) return false; 
    
    std::string line;
    while(std::getline(in, line)) {
        if(line.empty() || line[0]=='#') continue;
        std::stringstream ss(line);
        std::array<double,3> v;
        if(ss >> v[0] >> v[1] >> v[2]) 
            out.push_back(v);
    }
    return !out.empty();
}

// Detect lattice vectors from position data
static bool detect_lattice_vectors(
    const std::vector<vector2_t>& positions,
    vector2_t& a1, vector2_t& a2,
    int& Lx, int& Ly)
{
    if(positions.size() < 3) return false;

    // Collect unique coordinates
    std::vector<double> xs, ys;
    for(const auto& p : positions) {
        xs.push_back(p[0]);
        ys.push_back(p[1]);
    }
    std::sort(xs.begin(), xs.end());
    std::sort(ys.begin(), ys.end());
    
    // Remove duplicates
    auto last_x = std::unique(xs.begin(), xs.end(), 
        [](double a, double b) { return std::abs(a-b) < 1e-6; });
    xs.erase(last_x, xs.end());
    
    auto last_y = std::unique(ys.begin(), ys.end(), 
        [](double a, double b) { return std::abs(a-b) < 1e-6; });
    ys.erase(last_y, ys.end());

    Lx = xs.size();
    Ly = ys.size();
    
    if(Lx * Ly != positions.size()) {
        std::cerr << "Grid mismatch: " << Lx << "x" << Ly 
                  << " != " << positions.size() << "\n";
        return false;
    }

    // Compute lattice constants
    if(Lx > 1) {
        a1[0] = xs[1] - xs[0];
        a1[1] = 0.0;
    } else {
        a1 = {1.0, 0.0};  // Default if insufficient points
    }

    if(Ly > 1) {
        a2[0] = 0.0;
        a2[1] = ys[1] - ys[0];
    } else {
        a2 = {0.0, 1.0};  // Default if insufficient points
    }

    return true;
}

// Reorder positions/spins in row-major grid order
static void reorder_to_grid(
    const std::vector<vector2_t>& positions,
    const std::vector<Vectormath::Vector3>& spins,
    std::vector<vector2_t>& sorted_pos,
    std::vector<Vectormath::Vector3>& sorted_spins)
{
    std::vector<size_t> indices(positions.size());
    std::iota(indices.begin(), indices.end(), 0);
    
    std::sort(indices.begin(), indices.end(), 
        [&](size_t i, size_t j) {
            if(std::abs(positions[i][1] - positions[j][1]) < 1e-6)
                return positions[i][0] < positions[j][0];
            return positions[i][1] < positions[j][1];
        });
    
    sorted_pos.clear();
    sorted_spins.clear();
    sorted_pos.reserve(positions.size());
    sorted_spins.reserve(spins.size());
    
    for(auto i : indices) {
        sorted_pos.push_back(positions[i]);
        sorted_spins.push_back(spins[i]);
    }
}

extern "C" int compute_layered_topo_charge(
    const char* pos_file,
    const char* spin_file)
{
    std::vector<std::array<double,3>> pos3d, spin3d;
    if(!load_txt(pos_file, pos3d) || !load_txt(spin_file, spin3d)) {
        std::cerr << "Failed to load input files\n";
        return 1;
    }
    
    if(pos3d.size() != spin3d.size()) {
        std::cerr << "Position/spin count mismatch: " 
                  << pos3d.size() << " vs " << spin3d.size() << "\n";
        return 1;
    }

    // Group atoms by layer with tolerance relative to system height
    double min_z = 1e100, max_z = -1e100;
    for(const auto& p : pos3d) {
        if(p[2] < min_z) min_z = p[2];
        if(p[2] > max_z) max_z = p[2];
    }
    double z_tol = 1e-3 * (max_z - min_z);
    if(z_tol < 1e-10) z_tol = 1e-5;

    std::map<long, std::vector<int>> layer_map;
    for(int i=0; i<pos3d.size(); ++i) {
        long z_bin = std::lround(pos3d[i][2] / z_tol);
        layer_map[z_bin].push_back(i);
    }

    for(auto& [z_bin, indices] : layer_map) {
        double z_val = z_bin * z_tol;
        int N = indices.size();
        
        // Extract layer data
        std::vector<vector2_t> pos2d;
        std::vector<Vectormath::Vector3> spins;
        for(int i : indices) {
            pos2d.push_back({pos3d[i][0], pos3d[i][1]});
            spins.push_back({spin3d[i][0], spin3d[i][1], spin3d[i][2]});
        }

        // Detect lattice structure
        vector2_t a1, a2;
        int Lx, Ly;
        if(!detect_lattice_vectors(pos2d, a1, a2, Lx, Ly)) {
            std::cerr << "Skipping layer at z=" << z_val << "\n";
            continue;
        }

        // Reorder to grid layout
        std::vector<vector2_t> sorted_pos;
        std::vector<Vectormath::Vector3> sorted_spins;
        reorder_to_grid(pos2d, spins, sorted_pos, sorted_spins);

        // Generate triangulation using all points in the layer
        auto tris = compute_delaunay_triangulation_2D(sorted_pos);

        // Compute topological charge
        double Q = TopologicalCharge_FiniteLayer(sorted_spins, sorted_pos, tris);
        std::cout << "Layer z=" << z_val << " -> Q=" << Q << "\n";
    }
    return 0;
}

