// topological_charge.cpp
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <map>
#include <cmath>
#include <string>

// CGAL triangulation with info on vertices
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

using Kernel       = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point2       = Kernel::Point_2;
using Vb           = CGAL::Triangulation_vertex_base_with_info_2<std::size_t, Kernel>;
using Fb           = CGAL::Triangulation_face_base_2<Kernel>;
using Tds          = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using Triangulation= CGAL::Delaunay_triangulation_2<Kernel, Tds>;

struct Vec3 { double x, y, z;
    Vec3 operator-(Vec3 const& o) const { return {x-o.x, y-o.y, z-o.z}; }
    Vec3 cross(Vec3 const& o) const {
        return { y*o.z - z*o.y,
                 z*o.x - x*o.z,
                 x*o.y - y*o.x };
    }
    double dot(Vec3 const& o) const { return x*o.x + y*o.y + z*o.z; }
    Vec3 normalize() const {
        double n = std::sqrt(x*x + y*y + z*z);
        return {x/n, y/n, z/n};
    }
};

struct Triangle { int i0, i1, i2; };

// ——— Spirit‑style Delaunay triangulation over one unit cell ———
std::vector<Triangle> generate_spirit_triangles(
    std::vector<Vec3> const& positions,   // base‑cell atom positions
    int n_basis_atoms,
    int Lx, int Ly,
    bool pbc_x, bool pbc_y)
{
    // 1) Build base cell points: real atoms + three stretched corners
    double stretch = 0.1;
    Vec3 ta{1,0,0}, tb{0,1,0};
    int Nbase = n_basis_atoms + 3;
    std::vector<Vec3> base_pts(Nbase);
    for(int i=0; i<n_basis_atoms; ++i)
        base_pts[i] = {positions[i].x, positions[i].y, 0.0};

    // corner at a + b
    Vec3 bb = {ta.x + tb.x, ta.y + tb.y, 0.0};
    base_pts[n_basis_atoms]   = {positions[0].x + bb.x*stretch,
                                 positions[0].y + bb.y*stretch, 0.0};
    // corner at b
    base_pts[n_basis_atoms+1] = {positions[0].x + tb.x*(1-stretch),
                                 positions[0].y + tb.y*(1-stretch), 0.0};
    // corner at a
    base_pts[n_basis_atoms+2] = {positions[0].x + ta.x*(1-stretch),
                                 positions[0].y + ta.y*(1-stretch), 0.0};

    // 2) Triangulate base cell with CGAL, storing each vertex’s base index
    std::vector<std::pair<Point2,std::size_t>> pts;
    pts.reserve(Nbase);
    for(size_t i=0; i<base_pts.size(); ++i)
        pts.emplace_back(Point2(base_pts[i].x, base_pts[i].y), i);

    Triangulation dt;
    dt.insert(pts.begin(), pts.end());

    // 3) Extract base triangles
    std::vector<std::array<int,3>> base_tri;
    for(auto fit = dt.finite_faces_begin(); fit != dt.finite_faces_end(); ++fit) {
        base_tri.push_back({{
            int(fit->vertex(0)->info()),
            int(fit->vertex(1)->info()),
            int(fit->vertex(2)->info())
        }});
    }

    // 4) Translate each base-triangle across the Lx×Ly grid with PBC logic
    std::vector<Triangle> all_tri;
    all_tri.reserve(base_tri.size() * Lx * Ly);
    for(int b=0; b<Ly; ++b) {
        for(int a=0; a<Lx; ++a) {
            bool a_ok = (a+1 < Lx) || pbc_x;
            bool b_ok = (b+1 < Ly) || pbc_y;
            for(auto const& tri : base_tri) {
                std::array<int,3> idx;
                bool valid = true;
                for(int i=0; i<3; ++i) {
                    int vid = tri[i];
                    if(vid < n_basis_atoms) {
                        // same cell
                        idx[i] = vid + a*n_basis_atoms + b*n_basis_atoms*Lx;
                    }
                    else if(vid == n_basis_atoms && a_ok && b_ok) {
                        // translate by (a+1,b+1)
                        idx[i] = ((a+1)%Lx)*n_basis_atoms
                               + ((b+1)%Ly)*n_basis_atoms*Lx;
                    }
                    else if(vid == n_basis_atoms+1 && b_ok) {
                        // translate by (0,1)
                        idx[i] = a*n_basis_atoms
                               + ((b+1)%Ly)*n_basis_atoms*Lx;
                    }
                    else if(vid == n_basis_atoms+2 && a_ok) {
                        // translate by (1,0)
                        idx[i] = ((a+1)%Lx)*n_basis_atoms
                               + b*n_basis_atoms*Lx;
                    }
                    else {
                        valid = false;
                        break;
                    }
                }
                if(valid) all_tri.push_back({idx[0], idx[1], idx[2]});
            }
        }
    }

    return all_tri;
}

// ——— Berg–Lüscher solid angle ———
double solid_angle(Vec3 const& S1, Vec3 const& S2, Vec3 const& S3) {
    Vec3 a = S1.normalize();
    Vec3 b = S2.normalize();
    Vec3 c = S3.normalize();
    double triple = a.dot(b.cross(c));
    double denom  = 1.0 + a.dot(b) + b.dot(c) + c.dot(a);
    return 2.0 * std::atan2(triple, denom);
}

// ——— File loaders ———
bool load_spins(const std::string& fn, std::vector<Vec3>& spins) {
    std::ifstream in(fn);
    if(!in) return false;
    std::string line;
    while(std::getline(in, line)) {
        if(line.empty()||line[0]=='#') continue;
        std::stringstream ss(line);
        double sx, sy, sz;
        if(ss >> sx >> sy >> sz)
            spins.push_back({sx, sy, sz});
    }
    return true;
}

bool load_positions(const std::string& fn, std::vector<Vec3>& pos) {
    std::ifstream in(fn);
    if(!in) return false;
    std::string line;
    while(std::getline(in, line)) {
        if(line.empty()||line[0]=='#') continue;
        std::stringstream ss(line);
        double x,y,z;
        if(ss >> x >> y >> z)
            pos.push_back({x,y,z});
    }
    return true;
}

// ——— Slice into z‑layers ———
std::map<double,std::vector<std::pair<Vec3,Vec3>>> split_by_layers(
    std::vector<Vec3> const& positions,
    std::vector<Vec3> const& spins,
    double z_tol = 1e-6)
{
    std::map<double,std::vector<std::pair<Vec3,Vec3>>> layers;
    for(size_t i=0; i<positions.size(); ++i) {
        double z_key = std::round(positions[i].z / z_tol) * z_tol;
        layers[z_key].emplace_back(positions[i], spins[i]);
    }
    return layers;
}

// ——— Compute Q for one layer with Spirit‑style triangles ———
double compute_topological_charge_2D(
    std::vector<Vec3> const& pos2D,
    std::vector<Vec3> const& spins2D,
    int n_basis_atoms,
    int Lx, int Ly,
    bool pbc_x, bool pbc_y)
{
    auto triangles = generate_spirit_triangles(pos2D, n_basis_atoms, Lx, Ly, pbc_x, pbc_y);

     std::cout << "Generated triangles (" << triangles.size() << "):\n";
    for (size_t t = 0; t < triangles.size(); ++t) {
        const auto &tri = triangles[t];
        std::cout << "  [" << tri.i0 << ", " << tri.i1 << ", " << tri.i2 << "]\n";
    }



    double Q = 0.0, four_pi = 4.0 * M_PI;
    for(auto const& tri : triangles) {
        double ω = solid_angle(spins2D[tri.i0], spins2D[tri.i1], spins2D[tri.i2]);
        Q += ω / four_pi;
    }
    return Q;
}

// ——— C entry point for Python ctypes ———
extern "C"
void compute_layered_topo_charge(const char* pos_file, const char* spin_file) {
    std::vector<Vec3> positions, spins;
    if(!load_positions(pos_file, positions)) { std::cerr<<"Failed to load positions\n"; return; }
    if(!load_spins   (spin_file,   spins))    { std::cerr<<"Failed to load spins\n";    return; }
    if(positions.size()!=spins.size()) {
        std::cerr<<"Positions/spins size mismatch\n"; return;
    }

    auto layers = split_by_layers(positions, spins);
    for(auto const& [z, layer] : layers) {
        // Separate into 2D arrays
        std::vector<Vec3> pos2D, spin2D;
        for(auto const& [p, s] : layer) {
            pos2D.push_back(p);
            spin2D.push_back(s);
        }

        // For now assume 1-atom basis on a square Lx×Ly
        int total = int(pos2D.size());
        int L = int(std::round(std::sqrt(total)));
        if(L*L!=total) {
            std::cerr<<"Warning: non-square layer size="<< total <<"\n";
        }
        double Q = compute_topological_charge_2D(pos2D, spin2D,
            /*basis=*/1, L, L, /*pbc_x=*/false, /*pbc_y=*/false);
        std::cout<<"Layer z="<<z<<" → Q="<<Q<<"\n";
    }
}

// ——— Compile with:
//   g++ -O3 -std=c++17 -fPIC topological_charge.cpp -shared -o libtopo.so -lgmp -lmpfr

