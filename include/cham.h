#ifndef _CHAM_
#define _CHAM_

#include <cstring>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <fstream>
#include <set>
#include <algorithm>
#include <iomanip>
#include <list>
#include "outputMesh.h"

#pragma intrinsic(fabs)

// --- // ----------------------------------------------------- // --- //
// --- // ----------------------------------------------------- // --- //
//     //                    PLC CHAMFERING                     //     //
// --- // ----------------------------------------------------- // --- //
// --- // ----------------------------------------------------- // --- //

// #define PLCC_VERBOSE_DEBUG
  
//  "flat" edges will be ignored while aucte angles are searched.
//  "acute" if incident face form an acute dihedral angle.
//  "junk" edge have to be deleted during chamfering.
typedef enum class CHAMedge_t{
	undet, flat, acute, junk                               
} CHAMedge_type;

#define NO_BRIDGE UINT32_MAX

/// <summary>
/// CHAMedge
/// An edge or sub-edge chamfered PLC (see PLCc class). 
/// It is a segment.
/// ep[2] are the two endpoints of the edge.
/// In case the edge is a sub-edge (at least one endpint is a Steiner point)
/// oep[2] are the endpoints of the "older" parent-edge; otherwise oep[2]=ep[2].
/// <summary>
class CHAMedge{
public:
    uint32_t ep[2];                 // Endpoints (vertices-inds wrt PLCc 
                                    // vertices vector)
    uint32_t oep[2];                // "older" parent-edge endpoints 
                                    // (vertices-inds wrt PLCc vertices vector)
    std::vector<uint32_t> inc_face; // Incident CHAMfaces
    CHAMedge_type type;             // See CHAMedge_type note above.

    uint32_t loc_face_bridge_id; // If it is a "bridge-edge", contains the index 
                                 // of the related chamfered vertex; this info
                                 // is local wrt a face, i.e. bridges of  
                                 // different faces may have the same id.

    inline CHAMedge() : loc_face_bridge_id(NO_BRIDGE) {}
    inline CHAMedge(CHAMedge_type t) : type(t), loc_face_bridge_id(NO_BRIDGE) {}

    // copy constructor 
    inline CHAMedge(const CHAMedge& e) : 
        ep{e.ep[0], e.ep[1]}, oep{e.oep[0], e.oep[1]}, inc_face(e.inc_face), 
        type(e.type), loc_face_bridge_id(e.loc_face_bridge_id) {}

    // constructor for sub-edges
    inline CHAMedge(const uint32_t e0, const uint32_t e1,
                    const uint32_t oe0, const uint32_t oe1,
                    const std::vector<uint32_t>& ifaces, 
                    const CHAMedge_type t) : ep{e0, e1}, oep{oe0, oe1}, 
                                            inc_face(ifaces), type(t), 
                                            loc_face_bridge_id(NO_BRIDGE) {}

    // constructor for new edges
    inline CHAMedge(const uint32_t e0, const uint32_t e1, const uint32_t f) : 
                ep{e0, e1}, oep{e0, e1}, type(CHAMedge_t::undet), 
                loc_face_bridge_id(NO_BRIDGE) { 
        inc_face.push_back(f); 
    }

    inline void init_bridge_edge(const uint32_t e0, const uint32_t e1, 
                                 const uint32_t bridge_v, const uint32_t f){ 
        oep[0] = ep[0] = e0;  
        oep[1] = ep[1] = e1; 
        type = CHAMedge_t::undet; 
        loc_face_bridge_id = bridge_v; 
        inc_face.push_back(f); 
    }

    inline bool isFlat() const { return type == CHAMedge_t::flat; }
    inline bool isAcute() const { return type == CHAMedge_t::acute; }
    inline bool isJunk() const { return type == CHAMedge_t::junk; }
    inline bool isIsolated() const { return (inc_face.empty()); }
    inline bool isSubEdge() const { return (ep[0]!=oep[0] || ep[1]!=oep[1]); }
    inline bool isBridgeEdge() const { return loc_face_bridge_id != NO_BRIDGE; }

    inline bool hasVertex(uint32_t v) const { 
        return (ep[0] == v || ep[1] == v); 
    }
    inline bool hasOriginalVertex(uint32_t v) const { 
        return (oep[0] == v && oep[1] == v); 
    }
    bool coincident(const CHAMedge& e) const { 
        return (hasVertex(e.ep[0]) && hasVertex(e.ep[1])); 
    }
    bool has_commonVertex(const CHAMedge& e) const { 
        return commonVertex(e) != UINT32_MAX; 
    }

    uint32_t relatedOriginalVertex(uint32_t v) const {
        if (ep[0] == v) return oep[0];
        assert(ep[1] == v);
        return oep[1];
    }

    uint32_t oppositeVertex(uint32_t v) const {
        if (ep[0] == v) return ep[1];
        assert(ep[1] == v);
        return ep[0];
    }

    uint32_t oppositeOriginalVertex(uint32_t v) const {
        if (oep[0] == v) return oep[1];
        assert(oep[1] == v);
        return oep[0];
    }
    
    uint32_t commonVertex(const CHAMedge& e) const {
        if( e.hasVertex(ep[0]) ) return ep[0];
        if( e.hasVertex(ep[1]) ) return ep[1];
        return UINT32_MAX;
    }

    // Assumes that 'e' shares an endpoint with *this;
    // returs the endpoint of *tins not shared with 'e'
    uint32_t notCommonVertex(const CHAMedge& e) const {
        assert( has_commonVertex(e) );
        if( e.hasVertex(ep[0]) ) return ep[1];
                                 return ep[0];
    }

    // Returns TRUE if has the same loc_face_bridge_id
    bool same_bridge(const CHAMedge& e) const { 
        return loc_face_bridge_id == e.loc_face_bridge_id; 
    }

    // intialization 
    inline void set_ep(uint32_t e0, uint32_t e1){ ep[0] = e0; ep[1] = e1; }
    inline void set_ep_and_reset_oep(uint32_t e0, uint32_t e1){ 
        ep[0] = e0; ep[1] = e1; 
        reset_oep(); 
    }
    inline void reset_oep(){ oep[0] = ep[0]; oep[1] = ep[1]; }
    inline void reset_type(){ type = CHAMedge_t::undet; }

    // preEdge is an input triangle half-edge
    inline void fill_preEdge(uint32_t v0, uint32_t v1, uint32_t fi) {
        if (v0 > v1) { oep[0] = ep[0] = v1; oep[1] = ep[1] = v0; }
        else { oep[0] = ep[0] = v0; oep[1] = ep[1] = v1; }
        inc_face.push_back(fi);
    } 

    // modify icident faces
    // NOTE: 
    // removeIncidentFace and replaceIncidentFace do NOT update connectivity.

    void removeIncidentFace(uint32_t f) {
        for(uint32_t i=0; i<inc_face.size(); i++) if(inc_face[i]==f) { 
            inc_face.erase(inc_face.begin() + i); 
            return; 
        }
    }

    void replaceIncidentFace(uint32_t f_old, uint32_t f_new){
        for(uint32_t i=0; i<inc_face.size(); i++) if(inc_face[i]==f_old) { 
            inc_face[i]=f_new;            
            return; 
        }
    }

    void isolate(){ inc_face.clear(); }

    // Static functions to be used as predicates in std algorithms
    static inline bool isIsolatedPtr(const CHAMedge& e) {return e.isIsolated();}

    // Works as expected only if each edge has ep[0] < ep[1]
    static inline bool vertexSortFunc(const CHAMedge& e1, const CHAMedge& e2) {
        if (e1.ep[0] == e2.ep[0]) return (e1.ep[1] < e2.ep[1]);
        else return (e1.ep[0] < e2.ep[0]);
    }

    // DEBUG
    inline void print_parent_edge() const { 
        std::cout << "<" << oep[0] << "," << oep[1] << ">\n";
    }
};

inline std::ostream& operator<<(std::ostream& os, const CHAMedge& e) {
    os << "<" << e.ep[0] << "," << e.ep[1] << "> ";
    return os;
}


typedef std::vector<uint32_t>::iterator u32vect_iter;
#define EMPTY_PLACE UINT32_MAX

/// <summary>
/// CHAMface
/// This is a flat face of the chmfered PLC (see PLCc class)
/// It is an input triangle, or a portion of it (when chamfered)
/// When chamfering has been performed it has no acute angles.
/// </summary>
class CHAMface {
public:
    uint32_t triangle; // Index of the related input triangle
    std::vector<uint32_t> bounding_edges; // (Ordered) Set of bounding edges 
                                          // (represented as CHAMedge indices)

    CHAMface(){}
    CHAMface(uint32_t tri) : triangle(tri) {}

    inline bool isEmpty() const { return bounding_edges.empty(); }

    // Boundary navigation: works supposing that consecutive edges are stored 
    // consecutively in bounding_edges vector. 
    
    // index based
    void advance_on_bnd(size_t& it) const { 
        if( (++it) == bounding_edges.size() ) it = 0; 
    } 
    void advance_on_pierced_bnd(size_t& it) const { 
        size_t it0 = it;
        do{
            advance_on_bnd(it); 
            assert(it != it0);
        } while(bounding_edges[it]==EMPTY_PLACE); 
    } 
    void reverse_on_bnd(size_t& it) const { 
        if(it == 0) it = bounding_edges.size(); 
        it--; 
    }
    void reverse_on_pierced_bnd(size_t& it) const { 
        size_t it0 = it;
        do{
            reverse_on_bnd(it); 
            assert(it != it0);
        } while(bounding_edges[it]==EMPTY_PLACE); 
    } 
    void make_first(uint32_t ei){
        u32vect_iter it = std::find( bounding_edges.begin(), 
                                        bounding_edges.end(), ei);
        if(it!=bounding_edges.end()) make_first(it);
    }
    void make_last(uint32_t ei){
        u32vect_iter it = std::find(bounding_edges.begin(), 
                                        bounding_edges.end(), ei);
        if(it!=bounding_edges.end()) make_last(it);
    }

    // iterators based
    void advance_on_bnd(u32vect_iter& it) { 
        if( (++it) == bounding_edges.end() ) it = bounding_edges.begin(); 
    } 
    void reverse_on_bnd(u32vect_iter& it) { 
        if(it == bounding_edges.begin()) it = bounding_edges.end(); 
        it--; 
    }
    void make_first(u32vect_iter& it){ 
        std::rotate(bounding_edges.begin(), it, bounding_edges.end() ); 
    }
    void make_last(u32vect_iter& it){
        // To make the edge e (pointed by it) the last of bounding_edges,
        // rotate bounding_edges such that e+1 (the succesive on boundary of e) 
        // occupy the first position.
        advance_on_bnd(it);
        std::rotate( bounding_edges.begin(), it, bounding_edges.end() );
    }

    // boundary modification (NOTE: does not update connectivity).
    void replaceEdge_11(uint32_t old_e, uint32_t new_e){
       std::replace(bounding_edges.begin(), bounding_edges.end(), old_e, new_e);
    }

    // clean bounding_edges
    void remove_empty_places(){
        bounding_edges.erase(
                std::remove_if( bounding_edges.begin(), bounding_edges.end(), 
                                [](uint32_t x){ return x == EMPTY_PLACE; } ), 
            bounding_edges.end() );
    }

    // Static functions to be used as predicates in std algorithms
    static inline bool isEmptyPtr(const CHAMface& f) { return f.isEmpty(); }

    // only DEBUG functions
    inline bool hasBndEdge(uint32_t e_ind) const { 
        return (
            std::find( bounding_edges.begin(), bounding_edges.end(), e_ind) != 
            bounding_edges.end() ); 
    } 

    void print() const {
        std::cout<<"( "<< bounding_edges.size() <<" bound-edges )\n";
        for(size_t i = 0; i < bounding_edges.size(); i++){  
            std::cout<<"edge["<< bounding_edges[i] <<"]\n"; 
        }
        std::cout<<"\n";
    }
};

#define INVALID_BPT UINT32_MAX
// Chamfered PLC (removes angles < pi/2, from the input PLC)
class PLCc{
private:
    bool verbose;
    bool def_interior;
    bool manifold;
    bool safe_mode;

    const inputPLC& plc; // Input PLC enriched with stainer points + CDT

    double lfs;

public:
    const double epsilon;

    // input vertex -> incident triangle relaion;
    // computed in initalize() and used in search_acute_angles())
    std::vector< std::vector<uint32_t> > input_vt; 

    const size_t n_in_vrts;
    std::vector<pointType*> vertices; // Vertices of the chamfered PLC;
                                      // plc vertices (input vertices) are
                                    // stored in the first n_in_vrts positions.
    std::vector<CHAMedge> edges; // edges of the chamfered PLC
    std::vector<CHAMface> faces; // Faces of the chamfered PLC

    std::vector<double> vrt_ch_dist; // Contains chamfering (cut) distance for 
                                     // acute vertices, -1 otherwise.

    std::vector<uint32_t> ref_exp3D_vrt; // Has many entries as vertices vector; 
                                         // - for implicit vertices is the 
                                         // index of "closest" explicit3D point,
                                         // - for explicit points UINT32_MAX 

    // marker
    std::vector<uint32_t> mark_vrts;
    std::vector<uint32_t> mark_edges;
    std::vector<uint32_t> mark_faces;

    PLCc(const inputPLC& _plc, const double _epsilon, 
         bool _safe, bool simplify, bool _verbose) : plc(_plc), 
                                                     epsilon(_epsilon), 
                                                n_in_vrts(_plc.numVertices()), 
                                                     def_interior(true), 
                                                     manifold(true), 
                                                     safe_mode(_safe), 
                                                     verbose(_verbose),
                                                     lfs(DBL_MAX) {
        
        #ifdef PLCC_VERBOSE_DEBUG
        verbose = true;
        #endif

        if(verbose) std::cout<<"\nCHAMFERING:\n";
        initialize(); 
        search_acute_angles();    
        chamfering();
        if(simplify) chamfered_plc_simplification();
        if(verbose) std::cout<<"chamfering COMPLETED\n";
        assert( check_acuteness() );
    };

    // Access private fields
    inline bool input_plc_defines_interior() const { return def_interior; }
    inline bool input_plc_is_manifold() const { return manifold; }

    // 
    inline bool isSteinerVertex(uint32_t v) const { return v >= n_in_vrts; }
    inline bool is_acute_vrt(uint32_t vi) const { return vrt_ch_dist[vi] > 0.0;}
    inline bool is_acute_edge(uint32_t ei) const { return edges[ei].isAcute(); }
    inline bool has_acute_endpoint(uint32_t ei) const {
        return is_acute_vrt(edges[ei].ep[0]) || is_acute_vrt(edges[ei].ep[1]);
    }
    inline bool has_acute_edge(uint32_t fi) const {
        for(uint32_t be : faces[fi].bounding_edges ) if( be != EMPTY_PLACE ) {
            if( edges[ be ].isAcute() ) return true;
        }
        return false;
    }
    inline bool is_flat_edge(uint32_t ei) const { return edges[ei].isFlat(); }

    void get_face_vertices(const CHAMface& f, std::vector<uint32_t>& fv) const ;

    // Ordering
    void swap_edges(uint32_t e1, uint32_t e2); // swap edges[e1] and edges[e2] 
                                               // and update connectivity
    uint32_t move_back_isolated_edges();

    // Removes all isolated edges from edges vector
    void cleanUp_edges(){
        uint32_t edges_new_size = move_back_isolated_edges();
        edges.erase( std::next( edges.begin(), edges_new_size ), edges.end());
        mark_edges.resize(edges_new_size);
    }

    // Add new vertex to vertices vector and update related vectors,
    // returns the position of the added vertex in the vertices vector
    uint32_t add_vertex(pointType* p, uint32_t exp3d_i){
        uint32_t pos = (uint32_t)vertices.size();
        vertices.push_back( p );
        mark_vrts.push_back( 0 );
        vrt_ch_dist.push_back( -1.0 );
        ref_exp3D_vrt.push_back(exp3d_i);
        return pos;
    }

    // initialization - an input triangulated PLC is loaded
    //                - vertices, edges and triangular faces are created

    void initialize();
    void mergePreEdges(); // Removes duplicated pre-edges
    void orient_initial_triface_bnd(CHAMface& f);
    uint32_t inTri_opp_vrt(const CHAMedge& e, const uint32_t ti) const;
    void inTri_opp_edge(const uint32_t v, const uint32_t ti, 
                        uint32_t& u1, uint32_t& u2) const;

    bool findIF_flat_edge(const CHAMedge& e) const {
        if( e.inc_face.size()!=2 ) return false;
        return vOrient3D( e.ep[0], e.ep[1], 
                          inTri_opp_vrt(e, e.inc_face[0]), 
                          inTri_opp_vrt(e, e.inc_face[1]) ) == 0;
    }

    // Search for acute angle between PLCc elements
    void search_acute_angles();
    bool findIF_acute_edge(const CHAMedge& e) const;
    bool findIF_acute_vrt(uint32_t vi, const std::vector<uint32_t>& vv_i, 
                          const std::vector<uint32_t>& vt_i) const;

    double closest_vv_dist(uint32_t vi, const std::vector<uint32_t>& vv_i)const;
    double closest_vt_dist(uint32_t vi, const std::vector<uint32_t>& vt_i)const;
    double get_vrt_ch_dist(uint32_t vi, const std::vector<uint32_t>& vv_i, 
                           const std::vector<uint32_t>& vt_i) const {
        auto a = closest_vv_dist(vi, vv_i);
        auto b = closest_vt_dist(vi, vt_i);
        auto c = min(a, b);
        return c / 3;
    }
  
    // Chamfering
    void chamfering();
    void chamfering_vrts();
    void chamfer_edge_ep(const size_t ei, double d, const uint32_t ep_i);
    uint32_t new_vrt_on_segment(uint32_t v0, uint32_t v1, 
                                            double d, uint32_t d0);
    bool get_bpt_coeff(double& xi0, double& eta0, double& xi1, double& eta1, 
                        const vector3d& Ov, 
                        const vector3d& Ou0, const vector3d& Ou1, 
                        double t0, double t1,
                        const double zero_toll = 0.000000000000001);
    uint32_t new_vrts_in_inputTri(const uint32_t fi, const uint32_t vi, 
                                  const uint32_t u0, const uint32_t ou0, 
                                  const uint32_t u1, const uint32_t ou1);
    void chamfering_face(uint32_t fi);

    // Safe chamfering (optional)
    bool is_edge_normalizable(uint32_t ei) const {
        const pointType* ep0 = vertices[ edges[ei].ep[0] ];
        const pointType* ep1 = vertices[ edges[ei].ep[1] ];
        return (ep0->isLNC() && ep1->isBPT()) || (ep1->isLNC() && ep0->isBPT());
    }
    void normalize_face(uint32_t fi);

    // Simplification (Optional)
    typedef enum class bridge_simpl_t{
	    no, left, right, both                               
    } bridge_simpl_type;
    
    void chamfered_plc_simplification();
    uint32_t get_cons_edge_on_adj_face(uint32_t ei); // edges[ei] MUST have 
                                                     // only one incident face.
    bridge_simpl_t is_3bridge_simplifiable( uint32_t next_bel, uint32_t bel, 
                        uint32_t bec, uint32_t ber, uint32_t next_ber) const;
    void simplify_3bridge(uint32_t bel, uint32_t bec, uint32_t ber);
    void simplify_2bridge(uint32_t bes, uint32_t bec);
    // Final triangulation
    
    // Fills 'tri_fv' with a triangulation of the chamfered plc faces; each 
    // triangle is a triple of vertices inidces.
    void get_triangles(std::vector<uint32_t>& tri_fv) const; 
    void hear_clipping(uint32_t fi, 
                        std::vector<uint32_t>& out_tri_fv_list) const;

    
    // Given the triangulation of the chamfered plc 'out_tri', fills
    // 'compl_tri' with a triangulation of the complementar of the chamfered
    // plc wrt the input surface; each triangle is a triple of vertices inidces. 
    void get_complementar_tri(const std::vector<uint32_t>& out_tri, 
                                std::vector<uint32_t>& compl_tri);

    // Statistics
    double get_part_lfs() const {
        double shortest = DBL_MAX;
        for(const CHAMedge& e : edges) shortest = min(shortest, eEdgeSqLen(e)); 
        return min(sqrt(shortest), lfs);
    }

    // Predicates interfaces
    int vOrient3D(uint32_t v1, uint32_t v2, uint32_t v3, uint32_t v4) const {
        return -pointType::orient3D(*vertices[v1], *vertices[v2], 
                                    *vertices[v3], *vertices[v4]);
    }

    inline double vPtsSqDist(uint32_t v1, uint32_t v2) const{
        return ( vector3d(vertices[v1]) - vector3d(vertices[v2]) ).sq_length();
    }

    inline double eEdgeSqLen(const CHAMedge& e) const { 
        return vPtsSqDist( e.ep[0], e.ep[1] );
    }

    inline double eEdgeSqLen(const uint32_t ei) const { 
        return eEdgeSqLen(edges[ei]);
    }

    inline double vSqVrtDistSeg(uint32_t v, uint32_t u1, uint32_t u2) const {
        return vector3d(vertices[v]).sq_dist_segment(
            vector3d(vertices[ u1 ]), vector3d(vertices[ u2 ]));
    }

    inline double vSqVrtDistEdge(uint32_t v1, const CHAMedge& e) const {
        return vSqVrtDistSeg(v1, e.ep[0], e.ep[1]);
    }

    bool vCollinear(uint32_t v1, uint32_t v2, uint32_t v3) const {
        return !pointType::misaligned(*vertices[v1], 
                                      *vertices[v2], 
                                      *vertices[v3]);
    }

    bool vPointInInnerSegment(uint32_t u, uint32_t v1, uint32_t v2) const {
        return pointType::pointInInnerSegment(*vertices[u], 
                                              *vertices[v1], 
                                              *vertices[v2]);
    }

    bool vInnerSegmentCross(uint32_t u1, uint32_t u2, 
                            uint32_t v1, uint32_t v2) const {
        return pointType::innerSegmentsCross(*vertices[u1], *vertices[u2], 
                                             *vertices[v1], *vertices[v2]);
    }

    bool vInnerSegmentCrossesEdge(uint32_t u1, uint32_t u2, uint32_t ei) const {
        return vInnerSegmentCross(u1,u2,edges[ei].ep[0],edges[ei].ep[1]);
    }


    //Output

    void get_mesh_approximated_coords(std::vector<double>& approx_vrts) const {
        size_t num_vrts = vertices.size();
        approx_vrts.resize(num_vrts * 3);
        for(size_t i = 0; i < num_vrts; i++) {
            vertices[i]->getApproxXYZCoordinates( 
                approx_vrts[3*i], approx_vrts[3*i +1], approx_vrts[3*i +2] );
        }
    }

    void saveFaces() const {
        uint32_t num_vrts = (uint32_t)vertices.size();
        std::vector<double> out_coords;
        get_mesh_approximated_coords(out_coords);
        
        uint32_t nf = faces.size(); // number of cham faces.
        std::vector<uint32_t> nfv(nf); // number of vertices of each face.
        std::vector<uint32_t> fvind; // linearized vector of vertex indices of 
                                     // each face.
        for(size_t i = 0; i < faces.size(); i++) {
            std::vector<uint32_t> fv;  get_face_vertices(faces[i], fv);
            nfv[i] = (uint32_t)fv.size();
            fvind.insert(fvind.end(), fv.begin(), fv.end());
        }
        output_mesh cham_faces(out_coords.data(), num_vrts, fvind.data(), nf, 
                                "cut_plc_faces.off");
        cham_faces.save_polygon_mesh(nfv.data());
    }

    void saveFaces(std::vector<CHAMface>& save_faces, const char* title) const {
        uint32_t nv = (uint32_t)vertices.size();
        std::vector<double> out_coords;
        get_mesh_approximated_coords(out_coords);
        uint32_t nf = (uint32_t)save_faces.size();
        std::vector<uint32_t> fv;
        std::vector<uint32_t> nfv( save_faces.size() );
        uint32_t prec_cumulative_nfv = 0;
        for(size_t i = 0; i < save_faces.size(); i++) {
            get_face_vertices(save_faces[i], fv);
            nfv[i] = (uint32_t)fv.size() - prec_cumulative_nfv;
            prec_cumulative_nfv = (uint32_t)fv.size();
        }
        output_mesh cham_faces(out_coords.data(), nv, fv.data(), 1, title);
        cham_faces.save_minimal_polygon_mesh(nfv.data());
    }

    // save the triangles whose vertices indices are triple of consecutive 
    // integers of 'save' vector into a an .off file
    void saveTriangles(const std::vector<uint32_t>& tri, 
                                            const char* title) const {

        uint32_t nv = (uint32_t)vertices.size();
        std::vector<double> out_coords;
        get_mesh_approximated_coords(out_coords);

        assert( (tri.size() % 3) == 0 );

        uint32_t nt =  (uint32_t)tri.size()/3;
        output_mesh cham_tri_faces(out_coords.data(),nv,tri.data(),nt,title);
        cham_tri_faces.save_minimal_triangle_mesh();
    }

    void saveFace(const CHAMface& f, const char* title) const {

        std::vector<uint32_t> fv;
        get_face_vertices(f, fv);
        
        std::vector<uint32_t> nfv( { (uint32_t)fv.size() } );
        uint32_t nv = (uint32_t)vertices.size();
        std::vector<double> out_coords;
        get_mesh_approximated_coords(out_coords);
        output_mesh cham_faces(out_coords.data(), nv, fv.data(), 1, title);
        cham_faces.save_minimal_polygon_mesh(nfv.data());
    }

    void saveTriangle(uint32_t v1, uint32_t v2, uint32_t v3, 
                                                const char* title) const {
        double p1c[3], p2c[3], p3c[3];
        vertices[v1]->getApproxXYZCoordinates(p1c[0], p1c[1], p1c[2]);
        vertices[v2]->getApproxXYZCoordinates(p2c[0], p2c[1], p2c[2]);
        vertices[v3]->getApproxXYZCoordinates(p3c[0], p3c[1], p3c[2]);
        saveTriangle(p1c, p2c, p3c, title);
    }

    void saveTriangle(const double* p1c, const double* p2c, const double* p3c, 
                            const char* title) const {
        std::vector<double> out_coords( { p1c[0], p1c[1], p1c[2],
                                          p2c[0], p2c[1], p2c[2],
                                          p3c[0], p3c[1], p3c[2] });
        std::vector<uint32_t> tri( {0, 1, 2} );
        output_mesh out_tri(out_coords.data(),3, tri.data(),1, title);
        out_tri.save_minimal_triangle_mesh();
    }

    void saveEdgeIncFaces(const CHAMedge& e) const {
        std::vector<CHAMface> save_faces(e.inc_face.size());
        for (uint32_t fi : e.inc_face) save_faces[fi] = faces[ fi ];
        saveFaces(save_faces, "edge_inc_faces.off");
    }

    void save_rebuilded_input_after_chamfering(
            const std::vector<uint32_t>& tri_cham_faces, const char* title){
        std::vector<uint32_t> compl_tri;
        get_complementar_tri(tri_cham_faces, compl_tri);
		compl_tri.insert(compl_tri.end(), tri_cham_faces.begin(), tri_cham_faces.end());
		saveTriangles(compl_tri, title);
    }

    // DEBUG

    void print_vertex(uint32_t vi) const { 
        std::cout<<"vertex["<<vi<<"] = "<< *vertices[vi] <<"\n"; 
    }
    string vrt_type_toString(const pointType* v) const {
        if(v->isExplicit3D()) return "EXP3D";
        if(v->isLNC()) return "LNC";
        if(v->isBPT()) return "BPT";
        return "UNKNOWN";
    }
    string vrt_type_toString(uint32_t v) const { 
        return vrt_type_toString(vertices[v]);
    }
    void print_LNC_info(const pointType* v) const {
        const implicitPoint3D_LNC& l = v->toLNC();
        std::cout<<"( P = "<<l.P()<<", Q = "<<l.Q()<<", t = "<<l.T()<<" )\n";
    }
    void print_BPT_info(const pointType* v) const {
        const implicitPoint3D_BPT& b = v->toBPT();
        std::cout<<"( P = "<<b.P()<<", Q = "<<b.Q()<<", R = "<<b.R() 
                 <<", u = "<<b.U()<<", v = "<<b.V()<<" )\n";
    }
    void print_vertex_info(uint32_t vi) const { 
        const pointType* v = vertices[vi];
        print_vertex(vi);
        std::cout<<"TYPE: " << vrt_type_toString(vi) << " ";
        if(v->isLNC()) print_LNC_info(v);
        else if(v->isBPT()) print_BPT_info(v);
        else std::cout<<"\n";
    }

    inline void print_input_tri(uint32_t t) const {
        const uint32_t* tv = plc.triangle_vertices.data() + t*3;
        std::cout << *tv << " " << *(tv+1) << " " << *(tv+2) << "\n";
    }
    inline void print_input_vt(uint32_t vi) const { 
        std::cout<<"input vt["<<vi<<"]:\n";
        for(size_t i=0; i<input_vt[vi].size(); i++) 
            print_input_tri(input_vt[vi][i]);
        std::cout<<"\n";
    }

    void print_edge(size_t ei) const { 
        std::cout << "edge["<< ei <<"] = "<< edges[ei] <<"\n"; 
    }
    string edge_ep_types_toString(const CHAMedge& e) const {
        return string("(").append(
                    vrt_type_toString(e.ep[0])).append(",").append(
                    vrt_type_toString(e.ep[1])).append(") ");
    }
    inline void print_edge_and_epTypes(size_t ei) const { 
        std::cout <<"edge["<< ei <<"] = "<< edges[ei] <<" "
                  << edge_ep_types_toString(edges[ei]) << "\n";
    }
    void print_edge_and_inc_face(uint32_t ei) const { 
        std::cout<< "edge[" << ei << "] = "<< edges[ei] <<", inc_face: ";
        for(uint32_t fi : edges[ei].inc_face) std::cout<<fi<<" ";
        std::cout<<"\n"; 
    }

    void print_face(uint32_t fi) const {
        std::cout<<"face["<<fi<<"] : ";  faces[fi].print();
    }
    void print_face_edges(uint32_t fi) const {
        const std::vector<uint32_t>& fbnd = faces[fi].bounding_edges;
        std::cout<<"face["<<fi<<"] : ("<< fbnd.size() <<" bound-edges) \n";
        for(uint32_t ei : fbnd) 
            if(ei!=UINT32_MAX){ print_edge(ei); } else std::cout<<"JUNK edge\n";
        std::cout<<"\n";
    }

    bool check_vertex(uint32_t v) const {
        bool error = ( v >= vertices.size() || !vertices[v]->is3D());
        if(v >= vertices.size()) std::cout <<"vertex index "<< v <<" > vertices" 
                                           <<" size "<< vertices.size() <<"\n";
        if(!vertices[v]->is3D()) std::cout <<" vertex["<< v <<"] is not 3D\n"; 
        if(error) std::cout<<"[cham.h - vertex_has_sense()] ERROR dtected\n";
        return !error;
    }

    // During vertex chamfering (v0 is the chamfered vertex, the other 4 vi
    // are the ordered LNC-BPT-BPT-LNC bridge vertices)
    bool check_bridge(uint32_t v0, 
                    uint32_t v1, uint32_t v2, uint32_t v3, uint32_t v4) const {
        bool error = false;
        const uint32_t v[] = {v0, v1, v2, v3, v4};
        const pointType* po = vertices[v0];
        const pointType* p1 = vertices[v1];
        const pointType* p2 = vertices[v2];
        const pointType* p3 = vertices[v3];
        const pointType* p4 = vertices[v4];
        const bool v_match_type[] = { po->isExplicit3D(), p1->isLNC(),
                                        p2->isBPT(), p3->isBPT(), p4->isLNC() };
        const string type_str[] = {"EXPLICIT3D", "LNC", "BPT", "BPT", "LNC"};
        for(int i = 0; i < 5; i++) if( !v_match_type[i] ) {
            std::cout << "v" << i << " is not an " << type_str[i] << "\n"; 
            print_vertex_info(v[i]); std::cout<<"\n";
            error = true;
        }
        
        if(!error){
            const explicitPoint3D& expo = po->toExplicit3D();
            const pointType* p0 = (p1->toLNC().P().toExplicit3D() == expo) ? 
                                    &p1->toLNC().Q() : &p1->toLNC().P();
            const pointType* p5 = (p4->toLNC().P().toExplicit3D() == expo) ? 
                                    &p4->toLNC().Q() : &p4->toLNC().P();
            const pointType* p[] = {p0, p1, p2, p3, p4, p5};
            for(int i = 1; i < 5; i++) if(isAcuteAngle(p[i-1], p[i], p[i+1])) {
                std::cout<<"acute angle at v["<<i<<"]\n"; 
                error = true;
            }

            if(!error){
                if( pointType::innerSegmentsCross(*po, *p0, *p2, *p3) ||
                    pointType::innerSegmentsCross(*po, *p0, *p3, *p4) ||
                    pointType::innerSegmentsCross(*p1, *p2, *p3, *p4) ||
                    pointType::innerSegmentsCross(*po, *p5, *p1, *p2) ||
                    pointType::innerSegmentsCross(*po, *p5, *p2, *p3)   ) {
                    error = true;
                    std::cout<< "edge bridges overlap\n";
                }
            }
        }
        
        if(error) std::cout<<"[cham.h - check_bridge()] ERROR detected\n";
        return !error;
    }

    // during simplification
    void print_all_bridge_edges() const {
        std::cout<<"\nBRIDGE EDGES:\n";
        for(uint32_t ei=0; ei<edges.size(); ei++) if(edges[ei].isBridgeEdge()){ 
            const CHAMedge& e = edges[ei];
            std::cout <<"edge["<< ei <<"] = "<< edges[ei]
                    << edge_ep_types_toString(edges[ei]) << ", "
                    << "lfb_id = " << e.loc_face_bridge_id << ", " 
                    << "n_inc_faces = " << e.inc_face.size() <<"\n";
        }
        std::cout<<"\n";
    }

    void print_3bridge_details(uint32_t bel, uint32_t bec, uint32_t ber, 
                            uint32_t cons_left_edge, uint32_t cons_right_edge,
                                        bool print_inc_face_details) const {

        std::cout<<"\n";
        if(cons_left_edge==UINT32_MAX) std::cout<<"left_ext NOT FOUND!\n";
        else{ std::cout<<"left_ext: "; print_edge(cons_left_edge);}
        std::cout<<"left_bridge: "; print_edge(bel);  
        std::cout<<"cent_bridge: "; print_edge(bec);
        std::cout<<"right_bridge: "; print_edge(ber);
        if(cons_right_edge==UINT32_MAX) std::cout<<"right_ext NOT FOUND!\n";
        else{ std::cout<<"right_ext: "; print_edge(cons_right_edge);}
        if(!print_inc_face_details) return;
        assert(edges[ bec ].inc_face.size() == 1);
        const CHAMface& f = faces[ edges[ bec ].inc_face[0] ];
        for(size_t bei = 0; bei < f.bounding_edges.size(); bei++)
            print_edge_and_inc_face(f.bounding_edges[bei]);
    }

    void print_2bridge_details(uint32_t bec, uint32_t bes, 
                            uint32_t next_bes, bool print_inc_face_details) const {

        std::cout<<"\n";
        if(next_bes==UINT32_MAX) std::cout<<"left_ext NOT FOUND!\n";
        else{ std::cout<<"cons_side_ext: "; print_edge(next_bes);}
        std::cout<<"side_bridge: "; print_edge(bes);
        std::cout<<"cent_bridge: "; print_edge(bec);
        if(!print_inc_face_details) return;
        assert(edges[ bec ].inc_face.size() == 1);
        const CHAMface& f = faces[ edges[ bec ].inc_face[0] ];
        for(size_t bei = 0; bei < f.bounding_edges.size(); bei++)
            print_edge_and_inc_face(f.bounding_edges[bei]);
    }

    bool check_bridge_edges() const {
        for(const CHAMedge& e : edges) if(e.isBridgeEdge()) {
            if( e.inc_face.size() != 1) return false;
            const pointType* p0 = vertices[e.ep[0]];
            const pointType* p1 = vertices[e.ep[1]];
            if( !p0->isLNC() && !p0->isBPT() ) return false;
            if( !p1->isLNC() && !p1->isBPT() ) return false; 
            if( p0->isLNC() && p1->isLNC() ) return false; 
        }
        return true;
    }
    
    // Check that each acute edge has both acute endpoints: 
    // use it only before performing chamfering.
    bool acute_edges_have_acute_ep() const {
        bool error = false;
        for(uint32_t ei=0; ei<edges.size(); ei++) if(edges[ei].isAcute()){ 
            uint32_t n = 0;
            if(!is_acute_vrt(edges[ei].ep[0])) n++; 
            if(!is_acute_vrt(edges[ei].ep[1])) n++;
            if(n > 0){
                error = true;
                std::cout << "acute edge[" << ei << "] has "
                          << n << " not-acute endpoints\n";
            }
        }
        if(error) std::cout<<"[cham.h - acute_edges_have_acute_ep()] ERROR\n";
        return(!error);
    }

    bool check_subedge(uint32_t ei, uint32_t epi) const {
        const CHAMedge& e = edges[ei];
        uint32_t ep = e.ep[0], oep = e.oep[0], opp_ep = e.ep[1]; 
        if(epi==1){ ep = e.ep[1], oep = e.oep[1], opp_ep = e.ep[0]; }
        if(ep == oep) return true;
        bool error=false; 
        if(!vCollinear(oep, ep, opp_ep)){
            print_edge(ei);
            std::cout<<"sub-edge of "; e.print_parent_edge();
            std::cout<<"has "<<oep<<", "<<ep<<", "<<opp_ep<<" NOT collinear\n";
            error = true;
        }
        if(!vPointInInnerSegment(ep, oep, opp_ep)){
            print_edge(ei);
            std::cout<<"sub-edge of "; e.print_parent_edge();
            std::cout<<"has "<<ep<<" outside parent edge.\n";
            error = true;
        }
        return !error;
    }

    bool check_edge(uint32_t ei) const {
        const CHAMedge& e = edges[ei];
        bool error = false;
        if( ei >= edges.size()){
            error = true;
            std::cout<<"edge index "<<ei<<" > edges size "<<edges.size()<<"\n";
        }
        uint32_t e0 = e.ep[0], e1 = e.ep[1];
        if(!check_vertex(e0) || !check_vertex(e1)) error = true;
        if(e.isSubEdge() && 
            (!check_subedge(ei,0) || !check_subedge(ei,1)) ) error = true;
        if(error) std::cout << "[cham.h - check_edge()] ERROR detected\n";
        return !error;
    }

    bool vrt_inside_edge(uint32_t vi, uint32_t ei, bool verbose=false) const {
        bool inside = vPointInInnerSegment(vi, edges[ei].ep[0],edges[ei].ep[1]);
        if(inside && verbose) {
            std::cout << "vertex[" << vi << "] is inside "; print_edge(ei);
        }
        return inside;
    }

    bool edges_intersects(uint32_t ei, uint32_t ej, bool verbose=false) const { 
        uint32_t ei0 = edges[ei].ep[0], ei1 = edges[ei].ep[1];
        uint32_t ej0 = edges[ej].ep[0], ej1 = edges[ej].ep[1];
        if(vrt_inside_edge(ei0, ej, verbose)) return true;
        if(vrt_inside_edge(ei1, ej, verbose)) return true;
        if(vrt_inside_edge(ej0, ei, verbose)) return true;
        if(vrt_inside_edge(ej1, ei, verbose)) return true;
        bool intersects = vInnerSegmentCross(ei0, ei1, ej0, ej1);
        if(intersects && verbose) {
            std::cout<<"edge["<< ei <<"] and edge["<< ej <<"] intersects.\n";
            print_edge(ei); print_edge(ej); std::cout<<"\n";
        }
        return intersects;
    }

    bool face_boundary_overlaps(uint32_t fi, bool empty_places=false) const {
        const CHAMface& f = faces[fi];
        const std::vector<uint32_t>& fbe = f.bounding_edges;
        for(size_t i = 0; i < fbe.size()-1; i++) if(fbe[i] != EMPTY_PLACE) {
            for(size_t j = i+1; j < fbe.size(); j++) if(fbe[j] != EMPTY_PLACE) {
                if( !edges_intersects(fbe[i], fbe[j], true) ) continue;
                std::cout<<"[cam.h - face_boundary_overlaps()] ERROR: "
                         << "overlap detected on face "<< fi << "boundary\n";
                const uint32_t* tv = plc.triangle_vertices.data() + 
                                                    faces[fi].triangle * 3;
                saveTriangle(*tv, *(tv+1), *(tv+2), "overlapping_face_tri.off");
                return true;
            }       
        }
        return false;
    }

    bool check_face(uint32_t fi, bool empty_places=false) const {
        bool error = false;
        const CHAMface& f = faces[fi];
        const std::vector<uint32_t>& fbe = f.bounding_edges;

        // check tath each edge has at two neighbour (with which share an endpoint)
        for(size_t it=0; it < fbe.size(); ++it){

            if(empty_places && fbe[it] == EMPTY_PLACE) continue;

            size_t n_it = it;
            size_t p_it = it; 
            if(!empty_places) f.advance_on_bnd(n_it);
            else f.advance_on_pierced_bnd(n_it);
            if(!empty_places) f.reverse_on_bnd(p_it);
            else f.reverse_on_pierced_bnd(p_it);

            uint32_t ei = fbe[it];
            const CHAMedge& e = edges[ ei ];
            const CHAMedge& ne = edges[ fbe[n_it] ];
            const CHAMedge& pe = edges[ fbe[p_it] ];
            if(!e.has_commonVertex(ne) || !e.has_commonVertex(pe)){
                error = true;
                std::cout << "edge[" << ei << "] inchoerent bounding_edges\n";
            }
            if(e.coincident(ne) || e.coincident(pe)){
                error = true;
                std::cout << "edge[" << ei << "] bounding_edges duplication\n";
            }
        }

        if( face_boundary_overlaps(fi, empty_places) ) error = true;

        // check that each enpoint compares exactly two times on the bounday
        std::vector<uint32_t> loc_mark_vrts(vertices.size(), 0);
        for(uint32_t ei : f.bounding_edges) if(ei != EMPTY_PLACE) {
            loc_mark_vrts[edges[ei].ep[0]]++;
            loc_mark_vrts[edges[ei].ep[1]]++;
        }

        for(uint32_t ei : f.bounding_edges) if(ei != EMPTY_PLACE) 
            for(uint32_t k=0; k<2; k++) {
                uint32_t vi = edges[ei].ep[k];
                if(loc_mark_vrts[vi]!=2){
                    error = true;
                    std::cout << "endpoint " << vi << " of edge[" << ei 
                          << "] compares " << loc_mark_vrts[vi] << " times\n";
                }
            }

        if(error) {
            std::cout<<"[cham.h - check_face()] ERROR detected\n";
            print_face_edges(fi);
            return false;
        }
        return true;
    }

    void get_unchamferableFace_info(uint32_t fi) const {
        const uint32_t tri = faces[fi].triangle;
        const uint32_t* tv = plc.triangle_vertices.data() + tri*3;
        std::cout<<"chamfering failed for face["<< fi <<"] (input tri #"
                 << tri <<" <"<<*(tv)<<","<<*(tv+1)<<","<<*(tv+2)<<">)\n";
        print_face_edges(fi);
        saveTriangle(*tv, *(tv+1), *(tv+2), "unchamferable_tri.off");
    }

    bool check_EF_relation(uint32_t ei) const {
        const CHAMedge& e = edges[ei];
        for(uint32_t fi : e.inc_face) if(!faces[fi].hasBndEdge(ei)){ 
            std::cout<<"[cham.h - check_EF_relation()] ERROR "
                     <<"edge["<<ei<<"] has face["<<fi<<"] as inc_face, "
                     <<"but it does not appear on face boundary\n";
            return false; 
        }
        return true;
    }

    bool check_FE_relation(uint32_t fi) const{
        const CHAMface& f = faces[fi];
        for(uint32_t ei : f.bounding_edges){ 
            const CHAMedge& e = edges[ei];
            if(std::find(e.inc_face.begin(), e.inc_face.end(), fi) == e.inc_face.end()){ 
                std::cout<<"[cham.h - check_FE_relation()] ERROR face["<<fi<<"]"
                         <<" has edge["<<ei<<"] on its boundary, but it "
                        <<"does not appear as one of its incident faces\n";
                return false; 
            }
        }
        return true;
    }

    bool checkup() const {
        bool error = false;
        uint32_t ei, fi;
        for(ei=0; ei<edges.size();) if(!check_edge(ei++)) error=true;
        for(fi=0; fi<faces.size();) if(!check_face(fi++)) error=true;
        for(ei=0; ei<edges.size();) if(!check_EF_relation(ei++)) error = true;
        for(fi=0; fi<faces.size();) if(!check_FE_relation(fi++)) error = true;
        return !error;
    }

    // To check absence of acute angles after chamfering
    bool check_acuteness() const {

        bool error = false;
        const std::vector<pointType*>& v = vertices;

        // Check that each couple of faces incident at an edges do not form 
        // an acute dihedral angle
        for(uint32_t ei=0; ei<edges.size(); ei++){
            const CHAMedge& e = edges[ei];

            if(e.inc_face.size() == 1) continue;
            if(e.inc_face.size() > 4)
                std::cout<<"edges["<<ei<<"] is acute: it has "
                         << e.inc_face.size() <<" incident faces.\n";

            uint32_t ui, uj;
            const pointType* p0 = vertices[ e.oep[0] ];
            const pointType* p1 = vertices[ e.oep[1] ];
            assert(p0->isExplicit3D() && p1->isExplicit3D());
            const CHAMedge orig_e(e.oep[0], e.oep[1], UINT32_MAX);
            std::vector<uint32_t> in_tri;
            for(uint32_t fi : e.inc_face) in_tri.push_back(faces[fi].triangle);
            for(size_t i = 0; i < in_tri.size()-1; i++)
                for(size_t j = i+1; j < in_tri.size(); j++) {
                    ui = inTri_opp_vrt(orig_e, in_tri[i]);
                    uj = inTri_opp_vrt(orig_e, in_tri[j]);
                    if( isAcuteDihedral_exact(p0, p1, v[ui], v[uj]) ){
                        std::cout<<"edges["<<ei<<"] is acute.\n";
                        error = true;
                        // saveTriangle(ui,e.oep[0],o.oep[1], "acute_tri1.off");
                        // saveTriangle(uj,e.oep[0],o.oep[1], "acute_tri2.off");
                        // return false;
                    }
                }
            
        }

        // Check that each face has not acute angles at any two consecutive 
        // sides.
        for(uint32_t fi=0; fi<faces.size(); fi++){
            const CHAMface& f = faces[fi];
            const std::vector<uint32_t>& fbnd = f.bounding_edges;
            size_t it = 0;
            while(it < fbnd.size()){
                size_t nit = it; f.advance_on_bnd(nit);
                uint32_t vc, vl, vr;
                const CHAMedge& e1 = edges[ f.bounding_edges[it] ];
                const CHAMedge& e2 = edges[ f.bounding_edges[nit] ];
                vc = e1.commonVertex(e2);
                vl = e1.notCommonVertex(e2);
                vr = e2.notCommonVertex(e1);
                if( !e1.isFlat() && !e2.isFlat() &&
                    isAcuteAngle( v[vl], v[vc], v[vr]) ) {
                    std::cout<<"edge "<<e1<<" and edge "<<e2<<" form an "
                             <<"angle at face "<<fi<<" boundary.\n"; 
                    // saveFace(f, "acute_face.off");
                    // saveTriangle(vl, vc, vl, "acute_angle.off");
                    // return false;
                    error = true;
                }
                ++it;
            }
        }

        // Check that there are no acute vertices
        std::vector< std::vector<uint32_t> > ve_rel(vertices.size());
        for(size_t i = 0; i < edges.size(); i++) if(!edges[i].isFlat()) {
            ve_rel[ edges[i].ep[0] ].push_back(i);
            ve_rel[ edges[i].ep[1] ].push_back(i);
        }
        for(uint32_t vi=0; vi<vertices.size(); vi++) if(!ve_rel[vi].empty()) {
            assert(ve_rel[vi].size() > 1);
            for(size_t ei=0; ei<ve_rel[vi].size()-1; ei++)
                for(size_t ej=ei+1; ej<ve_rel[ei].size(); ej++) {
                    uint32_t ui = edges[ei].ep[0], uj = edges[ej].ep[0];
                    if(ui == vi) ui = edges[ei].ep[1];
                    if(uj == vi) uj = edges[ej].ep[1];
                    if(isAcuteAngle( v[ui], v[vi], v[uj]) ) {
                        std::cout<<"edge "<<ei<<" and edge "<<ej<<" form an "
                                 <<"angle at vertex "<<vi<<".\n"; 
                        // saveTriangle(ui, vi, vj, "acute_angle.off");
                        // return false;
                        error = true;
                    }
                }
        }

        if(error){
            std::cout<<"[cham.h - check_acuteness()] ERROR\n";
            return false;
        }
        return true;
    }

    //
    bool checkEdges_beforeJunkDeletion() const {
        for(const CHAMedge& e : edges)if(!e.isIsolated()){
            uint32_t e0 = e.ep[0], e1 = e.ep[1];
            if(e.isJunk()){ 
                if( !is_acute_vrt(e0) && !is_acute_vrt(e1) ) return false;
            }
            else if( is_acute_vrt(e0) || is_acute_vrt(e1) ) return false;
        }
        return true;
    }


    #ifdef PLCC_VERBOSE_DEBUG

    static inline bool sameLNC(const implicitPoint3D_LNC& l1, 
                                const implicitPoint3D_LNC& l2) {
        const explicitPoint3D& P1 = l1.P().toExplicit3D();
        const explicitPoint3D& Q1 = l1.Q().toExplicit3D();
        bigfloat t1 = l1.T();
        const explicitPoint3D& P2 = l2.P().toExplicit3D();
        const explicitPoint3D& Q2 = l2.Q().toExplicit3D();
        bigfloat t2 = l2.T();
        if(P1==P2 && Q1==Q2 && t1==t2) return true;
        if(P1==Q2 && Q1==P2 && t1==(1-t2)) return true;
        return false;
    }

    static inline bool sameBPT(const implicitPoint3D_BPT& b1, 
                                const implicitPoint3D_BPT& b2) {
        const explicitPoint3D& P1 = b1.P().toExplicit3D();
        const explicitPoint3D& Q1 = b1.Q().toExplicit3D();
        const explicitPoint3D& R1 = b1.R().toExplicit3D();
        bigfloat v1 = b1.V(), u1 = b1.U();
        const explicitPoint3D& P2 = b2.P().toExplicit3D();
        const explicitPoint3D& Q2 = b2.Q().toExplicit3D();
        const explicitPoint3D& R2 = b2.R().toExplicit3D();
        bigfloat v2 = b2.V(), u2 = b2.U();
        if(P1==P2 && Q1==Q2 && R1==R2 && v1==v2 && u1==u2) return true;
        if(P1==Q2 && Q1==P2 && R1==R2 && v1==u2 && u1==v2) return true;
        if(P1==R2 && Q1==P2 && R1==Q2 && v1==u2 && u1==(1-u2-v2)) return true;
        if(P1==P2 && Q1==R2 && R1==Q2 && v1==(1-u2-v2) && u1==u2) return true;
        if(P1==Q2 && Q1==R2 && R1==P2 && v1==(1-u2-v2) && u1==v2) return true;
        if(P1==R2 && Q1==Q2 && R1==P2 && v1==v2 && u1==(1-u2-v2)) return true;
        return false;
    }

    bool duplicated_vertices() const {
        bool duplicated = false;
        const std::vector<pointType*>& v = vertices;
        std::vector< std::pair<vector3d,uint32_t> > v_cpy;
        for(size_t i=0; i<vertices.size(); i++)
            v_cpy.push_back( std::pair<vector3d,uint32_t> (vector3d(v[i]),i) ); 
        
        std::sort(v_cpy.begin(), v_cpy.end(),
                    []( const std::pair<vector3d,uint32_t> &a, 
                        const std::pair<vector3d,uint32_t> &b ) { 
                            return ( a.first < b.first || 
                                (a.first == b.first && a.second < b.second) ); 
                        } );
        
        for(size_t i=1; i<v_cpy.size(); i++) {
            if( !(v_cpy[i-1].first == v_cpy[i].first) ) continue;
            
            const pointType* p1 = v[ v_cpy[i-1].second ];
            const pointType* p2 = v[ v_cpy[i].second ];
            if( p1->isExplicit3D() == p2->isExplicit3D() ||
               (    p1->isLNC() && p2->isLNC() && 
                    sameLNC(p1->toLNC(),p2->toLNC())   ) ||
               (    p1->isBPT() && p2->isBPT() && 
                    sameBPT(p1->toBPT(),p2->toBPT())   )    ){  
                duplicated = true;
                std::cout << "duplicated vertex[" << v_cpy[i-1].second
                          << "] = vertex[" << v_cpy[i].second << "]\n"; 
                print_vertex_info(v_cpy[i-1].second); 
                print_vertex_info(v_cpy[i].second);
            }
        }
        if(duplicated) std::cout<<"[cham.h - duplicated_vertices()] ERROR\n";
        return duplicated;
    }

    bool duplicated_edges() const{
        bool duplicated = false;
        std::vector< std::vector<uint32_t> > edge_cpy;
        for(uint32_t i=0; i<edges.size(); i++){ 
            uint32_t ep0=edges[i].ep[0], ep1=edges[i].ep[1];
            if(ep0 > ep1) std::swap(ep0, ep1);
            edge_cpy.push_back( std::vector<uint32_t> ({ep0, ep1, i}) ); 
        }
        std::sort(edge_cpy.begin(), edge_cpy.end(),
                  []( const std::vector<uint32_t> &a, 
                      const std::vector<uint32_t> &b ) { 
                        return (a[0] < b[0] || (a[0]==b[0] && a[1]<b[1])); } );
        for(size_t i=1; i<edge_cpy.size(); i++){ 
            if( edge_cpy[i-1][0] == edge_cpy[i][0] &&
                edge_cpy[i-1][1] == edge_cpy[i][1]    ){
                std::cout<<"duplicated edge["<< edge_cpy[i-1][2] 
                        <<"] = edge["<< edge_cpy[i][2] <<"]\n"; 
                duplicated = true;
            }
        }
        if(duplicated) std::cout<<"[cham.h - duplicated_edges()] ERROR\n";
        return duplicated;
    }

    inline void disp_howManyFlatEdges() const {
        uint32_t count_flat = 0;
        for(const CHAMedge& e : edges) if( e.isFlat() ) count_flat++;
        std::cout<<count_flat<<" flat edges.\n";
    }

    inline void disp_howManyAcuteEdges() const {
        uint32_t count_acute_e = 0;
        for(const CHAMedge& e : edges) if(e.isAcute()) count_acute_e++;
        std::cout<<count_acute_e<<" acute edges.\n";
    }

    inline void disp_howManyAcuteVertices() const {
        uint32_t count_acute_v = 0;
        for(uint32_t i=0; i<vertices.size(); i++) if(is_acute_vrt(i)) count_acute_v++;
        std::cout<<count_acute_v<<" acute vertices founded.\n";
    }

    inline void print_just_chamfered_edge(uint32_t ei, uint32_t ep_i) const {
        std::cout<<"after chamfering:\n"; 
        print_edge(ei); print_edge(edges.size()-1);
        std::cout<<"edge["<< ((ep_i==1) ? edges.size()-1 : ei) <<"] is junk\n";
    }

    #endif // PLCC_VERBOSE_DEBUG

    void report_error_and_exit(const char* msg) const {
        printf("%s", msg);
        exit(1);
    }

};

#include "cham.cpp"

#endif // _CHAM_
