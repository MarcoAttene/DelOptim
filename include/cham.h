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

#pragma intrinsic(fabs)

// --- // ----------------------------------------------------- // --- //
// --- // ----------------------------------------------------- // --- //
//     //                    PLC CHAMFERING                     //     //
// --- // ----------------------------------------------------- // --- //
// --- // ----------------------------------------------------- // --- //

//#define PLCC_DEBUG
//#define PLCC_VERBOSE_DEBUG
//#define PLCC_VERBOSE_DEBUG_LEV1

//#define FP_CHAMFERING

#ifdef TEST_CHAMFERING
// Exit values for testing purposes
typedef enum class EXIT_t{
	ok, open_input, expected_fail, vons_fail                  
} EXIT_type;
#endif

// NOTES:   
//  "flat" edges will be ignored.
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
/// oep[2] are the endpoints of the "older" parent-edge; otherwise oep[2] = ep[2].
/// <summary>
class CHAMedge{
public:
    uint32_t ep[2];                 // Endpoints (vertices-inds wrt PLCc vertices vector)
    uint32_t oep[2];                // "older" parent-edge endpoints (vertices-inds wrt PLCc vertices vector)
    std::vector<uint32_t> inc_face; // Incident CHAMfaces
    CHAMedge_type type;             // See CHAMedge_type note above.

    uint32_t loc_face_bridge_id; // If it is a "bridge-edge" contains the index of the related chamfered vertex
                                 // this info is local wrt a face, i.e. bridges of diffent faces may have the same id.

    inline CHAMedge() : loc_face_bridge_id(NO_BRIDGE) {}
    inline CHAMedge(CHAMedge_type t) : type(t), loc_face_bridge_id(NO_BRIDGE) {}

    // copy constructor 
    inline CHAMedge(const CHAMedge& e) : ep{e.ep[0], e.ep[1]}, oep{e.oep[0], e.oep[1]}, inc_face(e.inc_face), type(e.type), loc_face_bridge_id(e.loc_face_bridge_id) {}

    // constructor for sub-edges
    inline CHAMedge(const uint32_t e0, const uint32_t e1, const uint32_t oe0, const uint32_t oe1,
        const std::vector<uint32_t>& ifaces, const CHAMedge_type t) : ep{e0, e1}, oep{oe0, oe1}, inc_face(ifaces), type(t), loc_face_bridge_id(NO_BRIDGE) {}

    // constructor for new edges
    inline CHAMedge(const uint32_t e0, const uint32_t e1, const uint32_t f) : ep{e0, e1}, oep{e0, e1}, type(CHAMedge_t::undet), loc_face_bridge_id(NO_BRIDGE) { inc_face.push_back(f); }

    inline void init_bridge_edge(const uint32_t e0, const uint32_t e1, const uint32_t bridge_v, const uint32_t f){ 
        oep[0] = ep[0] = e0;  oep[1] = ep[1] = e1; 
        type = CHAMedge_t::undet; 
        loc_face_bridge_id = bridge_v; 
        inc_face.push_back(f); 
    }

    // Enquiry
    inline bool isFlat() const { return type == CHAMedge_t::flat; }
    inline bool isAcute() const { return type == CHAMedge_t::acute; }
    inline bool isJunk() const { return type == CHAMedge_t::junk; }
    inline bool isIsolated() const { return (inc_face.empty()); }

    inline bool hasVertex(uint32_t v) const { return (ep[0] == v || ep[1] == v); }
    inline bool hasOriginalVertex(uint32_t v) const { return (oep[0] == v && oep[1] == v); }
    bool coincident(const CHAMedge& e) const { return (hasVertex(e.ep[0]) && hasVertex(e.ep[1])); }
    bool has_commonVertex(const CHAMedge& e) const { return commonVertex(e)!=UINT32_MAX; }

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

    // assumes that e shares an endpoint with *this
    // returs the endpoint of not shared with e
    uint32_t notCommonVertex(const CHAMedge& e) const {
        assert( has_commonVertex(e) );
        if( e.hasVertex(ep[0]) ) return ep[1];
                                 return ep[0];
    }

    // intialization 

    inline void set_ep(uint32_t e0, uint32_t e1){ ep[0]=e0; ep[1]=e1; }
    inline void set_ep_and_reset_oep(uint32_t e0, uint32_t e1){ ep[0]=e0; ep[1]=e1; reset_oep(); }
    inline void reset_oep(){ oep[0]=ep[0]; oep[1]=ep[1]; }
    inline void reset_type(){ type = CHAMedge_t::undet; }

    // preEdge is an input triangle half-edge
    inline void fill_preEdge(const uint32_t v0, const uint32_t v1, const uint32_t fi) {
        if (v0 > v1) { oep[0] = ep[0] = v1; oep[1] = ep[1] = v0; }
        else { oep[0] = ep[0] = v0; oep[1] = ep[1] = v1; }
        inc_face.push_back(fi);
    } 

    // modify icident faces
    // NOTE: removeIncidentFace and replaceIncidentFace functions do not update connectivity

    void removeIncidentFace(uint32_t f) {
        for(uint32_t pos=0; pos<inc_face.size(); pos++) if(inc_face[pos]==f) { 
            inc_face.erase(inc_face.begin() + pos); 
            return; 
        }
    }

    void replaceIncidentFace(uint32_t f_old, uint32_t f_new){
        for(uint32_t pos=0; pos<inc_face.size(); pos++) if(inc_face[pos]==f_old) { 
            inc_face[pos]=f_new;            
            return; 
        }
    }

    void isolate(){ inc_face.clear(); }

    // Static functions to be used as predicates in std algorithms
    static inline bool isIsolatedPtr(const CHAMedge& e) { return e.isIsolated(); }

    // Works as expected only if each edge has ep[0] < ep[1]
    static inline bool vertexSortFunc(const CHAMedge& e1, const CHAMedge& e2) {
        if (e1.ep[0] == e2.ep[0]) return (e1.ep[1] < e2.ep[1]);
        else return (e1.ep[0] < e2.ep[0]);
    }
};


typedef std::vector<uint32_t>::iterator u32vect_iter;
#define EMPTY_PLACE UINT32_MAX

/// <summary>
/// CHAMface
/// This is a maximal flat face of the chmfered PLC (see PLCc class)
/// It is built out of edge-adjacent and coplanar parts of input triangles
/// It has no acute corners
/// It might be bounded by one or more loops of CHAMedges
/// </summary>
class CHAMface {
public:
    uint32_t triangle; // One original triangle coplanar and 2D-intersecting the face: 
                       // it carryes the the orientation wrt the input surface.
    std::vector<uint32_t> bounding_edges; // (Ordered) Set of bounding edges (CHAMedge indices)
    bool is_simply_connected = true;

    CHAMface(){}
    CHAMface(uint32_t tri) : triangle(tri) {}

    inline bool isEmpty() const { return bounding_edges.empty(); }

    inline bool hasBndEdge(uint32_t e_ind) const { return (std::find( bounding_edges.begin(), bounding_edges.end(), e_ind) != bounding_edges.end()); } // DEBUG

    // boundary navigation: 
    // works supposing that boundary edges forms a unique chain (ring) when the face is simply connected,
    // and that boundary edges of not-simply connected faces form chains that are stored contiguosly. 
    void advance_on_bnd(size_t& it) const { it++; if( it == bounding_edges.size() ) it=0; } 
    void reverse_on_bnd(size_t& it) const { if(it==0) it = bounding_edges.size(); it--; }
    void advance_on_bnd(u32vect_iter& it) { it++; if( it == bounding_edges.end() ) it = bounding_edges.begin(); } 
    void reverse_on_bnd(u32vect_iter& it) { if(it == bounding_edges.begin()) it = bounding_edges.end(); it--; }
    void make_first(u32vect_iter& it){ std::rotate(bounding_edges.begin(), it, bounding_edges.end() ); }
    void make_first(uint32_t ei);
    void make_last(u32vect_iter& it);
    void make_last(uint32_t ei);
    void replaceEdge_11(uint32_t old_e, uint32_t new_e); // NOTE: this function does not update connectivity.

    // Static functions to be used as predicates in std algorithms
    static inline bool isEmptyPtr(const CHAMface& f) { return f.isEmpty(); }
};

#define INVALID_BPT UINT32_MAX
// Chamfered PLC (remove from input PLC angles < pi/2)
class PLCc{
private:
    bool verbose;
    bool def_interior;
    bool manifold;

public:
    const double epsilon;

    const inputPLC& plc; // Input PLC enriched with stainer points + CDT
    std::vector<std::vector<uint32_t>> input_vt; // input vertex->incident_triangle relaion (computed in initalize() and used in search_acute_angles())

    const size_t n_in_vrts;
    std::vector<pointType*> vertices; // vertices of the chamfered PLC
                                        // first plcx.input_nv vertices (input vertices) 
                                        // have to be stored in the first n_in_vrts positions.
    std::vector<CHAMedge> edges; // edges of the chamfered PLC
    std::vector<CHAMface> faces; // Faces of the chamfered PLC

    std::vector<double> vrt_ch_dist;  // contains chamfering (cut) distance if v is an acute vertex, -1 otherwise  

    std::vector<uint32_t> ref_exp3D_vrt; // store the vrtex index of "closest" explicit3D point for implicit points, UINT32_MAX for explicit points

    // marker
    std::vector<uint32_t> mark_vrts;
    std::vector<uint32_t> mark_edges;
    std::vector<uint32_t> mark_faces;

    PLCc(const inputPLC& _plc, const double _epsilon, bool _verbose) : plc(_plc), epsilon(_epsilon), n_in_vrts(_plc.numVertices()), def_interior(true), manifold(true), verbose(_verbose) {
        
        #ifdef PLCC_VERBOSE_DEBUG
        verbose = true;
        #endif

        if(verbose) std::cout<<"\nCHAMFERING:\n";

        initialize(); 
        if(verbose) std::cout<<"[PLCc] - initialization COMPLETED\n";

        search_acute_angles();

        if(verbose) std::cout<<"[PLCc] - determination of acute vertices and edges COMPLETED\n";

        #ifdef PLCC_DEBUG
        assert( acute_edges_have_acute_ep() && checkup() );
        std::cout<<"[PLCc] - pre-chamfering debug COMPLETED\n";
        #endif
        
        chamfering();

        #ifdef PLCC_DEBUG
        assert( checkup() );
        std::cout<<"[PLCc] - post-chamfering debug COMPLETED\n";
        #endif
        
        if(verbose) std::cout<<"[PLCc] - chamfering COMPLETED\n";

        // if( !check_acuteness() ) exit(99); // DEBUG
    };

    inline bool isSteinerVertex(uint32_t v) const { return v >= n_in_vrts; }
    inline bool is_acute_vrt(const uint32_t vi) const { return vrt_ch_dist[vi] > 0.0; }
    inline bool is_acute_edge(const uint32_t ei) const { return edges[ei].isAcute(); }
    inline bool is_flat_edge(const uint32_t ei) const { return edges[ei].isFlat(); }

    void get_face_vertices(const CHAMface& f, std::vector<uint32_t>& fv) const ;

    // Ordering

    void swap_edges(uint32_t e1, uint32_t e2); // swap edges[e1] and edges[e2] and update connectivity
    uint32_t move_back_isolated_edges();

    // Removes all isolated edges from edges vector
    void cleanUp_edges(){
        uint32_t edges_new_size = move_back_isolated_edges();
        edges.erase( std::next( edges.begin(), edges_new_size ), edges.end());
        mark_edges.resize(edges_new_size);
    }

    // Add new vertex to vertices vector and update related vectors
    inline void add_vertex(pointType* p, uint32_t exp3d_i){
        vertices.push_back( p );
        mark_vrts.push_back( 0 );
        vrt_ch_dist.push_back( -1.0 );
        ref_exp3D_vrt.push_back(exp3d_i);
    }

    // void swap_faces(uint32_t f1, uint32_t f2); // swap faces[f1] and faces[f2] and update connectivity
    // uint32_t move_back_empty_faces();

    // initialization - an input triangulated PLC is loaded
    //                - vertices, edges and triangular faces are created
    //                - ajacent coplanar faces are merged: flat edges and empty faces (due to merging) are removed.

    void initialize();

    void mergePreEdges(); // Removes duplicated pre-edges

    void orient_initial_triface_bnd(CHAMface& f);
    uint32_t inTri_opp_vrt(const CHAMedge& e, const uint32_t ti) const;
    void inTri_opp_edge(const uint32_t v, const uint32_t ti, uint32_t& u1, uint32_t& u2) const;

    bool findIF_flat_edge(const CHAMedge& e) const {
        if( e.inc_face.size()!=2 ) return false;
        return vOrient3D(e.ep[0], e.ep[1], inTri_opp_vrt(e, e.inc_face[0]), inTri_opp_vrt(e, e.inc_face[1])) == 0;
    }

    bool input_plc_defines_interior() const { return def_interior; }
    bool input_plc_is_manifold() const { return manifold; }

    // Search for acute angle between PLCc elements
    
    void search_acute_angles();

    bool findIF_acute_edge(const CHAMedge& e) const;
    bool findIF_acute_vrt(const uint32_t vi, const std::vector<uint32_t>& vv_i, const std::vector<uint32_t>& vt_i) const;

    double closest_vv_dist(const uint32_t vi, const std::vector<uint32_t>& vv_i) const;
    double closest_vt_dist(const uint32_t vi, const std::vector<uint32_t>& vt_i) const;

    double get_vrt_ch_dist(uint32_t vi, const std::vector<uint32_t>& vv_i, const std::vector<uint32_t>& vt_i) const {
        auto a = closest_vv_dist(vi, vv_i);
        auto b = closest_vt_dist(vi, vt_i);
        auto c = min(a, b);
        return c / 3;
    }

    // Input triangle chamfering (Fallback solution)

    void inputTriangleChamfering(uint32_t fi);
  
    // Chamfering

    void chamfering();
    void chamfering_vrts();
    void chamfer_edge_ep(const size_t ei, double d, const uint32_t ep_i);
    uint32_t new_vrt_on_segment(uint32_t v0, uint32_t v1, const double d, const uint32_t d0);
    uint32_t new_vrts_in_inputTri(const uint32_t fi, const uint32_t vi, const uint32_t u0, const uint32_t ou0, const uint32_t u1, const uint32_t ou1);
    void remove_junk_edges_from_face(uint32_t fi);
    void get_edge_ch_dist(std::vector<double>& edge_ch_dist);
    void chamfer_acute_edge_from_inc_face(uint32_t ei, uint32_t fi, double d);
    uint32_t chamfering_face(uint32_t fi, bool safe);

    // Simplification (Post-Processing)

    void chamfered_plc_simplification();
    uint32_t get_cons_edge_on_adj_face(uint32_t ei); // edges[ei] MUST have only one incident face.
    
    // Final triangulation

    void hear_clipping(uint32_t fi, std::vector<uint32_t>& out_tri_fv_list) const;
    void triangulate_chamfered_plc(); // WORK IN PROGRESS

    // get a triengulation of the complementar of the chamfered faces wrt the input surface
    void get_complementar_tri(const std::vector<uint32_t>& out_tri, std::vector<uint32_t>& compl_tri);

    // Predicates interfaces

    // vertices index interface for geometric predicates
    int vOrient3D(uint32_t v1, uint32_t v2, uint32_t v3, uint32_t v4) const {
        return -pointType::orient3D(*vertices[v1], *vertices[v2], *vertices[v3], *vertices[v4]);
    }

    inline double vPtsSqDist(uint32_t v1, uint32_t v2) const{
        return ( vector3d( vertices[v1] ) - vector3d( vertices[v2] ) ).sq_length() ;
    }

    inline double eEdgeSqLen(const uint32_t ei) const { return vPtsSqDist( edges[ei].ep[0], edges[ei].ep[1] ); }

    inline double vSqVrtDistSeg(uint32_t v, uint32_t u1, uint32_t u2) const {
        return vector3d(vertices[v]).sq_dist_segment(vector3d(vertices[ u1 ]), vector3d(vertices[ u2 ]));
    }

    inline double vSqVrtDistEdge(uint32_t v1, const CHAMedge& e) const {
        return vSqVrtDistSeg(v1, e.ep[0], e.ep[1]);
    }

    bool vCollinear(uint32_t v1, uint32_t v2, uint32_t v3) const {
        return !pointType::misaligned(*vertices[v1], *vertices[v2], *vertices[v3]);
    }

    bool vPointInInnerSegment(uint32_t u, uint32_t v1, uint32_t v2) const {
        return pointType::pointInInnerSegment(*vertices[u], *vertices[v1], *vertices[v2]);
    }

    bool vInnerSegmentCross(uint32_t u1, uint32_t u2, uint32_t v1, uint32_t v2) const {
        return pointType::innerSegmentsCross(*vertices[u1], *vertices[u2], *vertices[v1], *vertices[v2]);
    }

    bool vInnerSegmentCrossesEdge(uint32_t u1, uint32_t u2, uint32_t ei) const {
        return vInnerSegmentCross(u1,u2,edges[ei].ep[0],edges[ei].ep[1]);
    }


    //Output

    void get_trivrts_1av(const std::vector<uint32_t>& fv, std::vector<uint32_t>& tri_fv) const;
    void get_trivrts_2av(const std::vector<uint32_t>& fv, std::vector<uint32_t>& tri_fv) const;
    void get_trivrts_2av1ae(const std::vector<uint32_t>& fv, std::vector<uint32_t>& tri_fv) const;
    void get_trivrts_3av(const std::vector<uint32_t>& fv, std::vector<uint32_t>& tri_fv) const;
    void get_trivrts_3av1ae(const std::vector<uint32_t>& fv, std::vector<uint32_t>& tri_fv) const;
    void get_trivrts_3av2ae(const std::vector<uint32_t>& fv, std::vector<uint32_t>& tri_fv) const;
    void get_trivrts_3av3ae(const std::vector<uint32_t>& fv, std::vector<uint32_t>& tri_fv) const;
    void get_triangles(std::vector<uint32_t>& tri_fv) const; // triangulate faces without taking in to account connectivity
    void get_triangles_naive(std::vector<uint32_t>& tri_fv) const; // DO NOT SUPPORT SIMPLICATION, triangulate faces without taking in to account connectivity

    // Export the data structure to optimizator
    void chamferPLC(inputPLC& _plc, double _epsilon, std::vector<genericPoint *>& _vertices, std::vector<std::vector<uint32_t>>& _faces){
        _plc.coordinates.assign(plc.coordinates.begin(), plc.coordinates.end());
        _plc.triangle_vertices.assign(plc.triangle_vertices.begin(), plc.triangle_vertices.end());
        _epsilon = epsilon;
        _vertices.assign(vertices.begin(), vertices.end());

        std::vector<uint32_t> out_tri;
        get_triangles(out_tri);
        _faces.resize(out_tri.size() /3);
        size_t out_fi = 0;
        for(size_t i=0; i<out_tri.size()/3; i++) 
            _faces[out_fi++].assign(out_tri.begin() + 3*i, out_tri.begin()+ 3*i +3);
    }

    // do not work if non-simpy connected faces are present
    bool saveFaces() const {

        for(const CHAMface& f : faces)if( !f.is_simply_connected ){  
            std::cout<<"PLCc - not simply connected face founded (off file not generated)\n";
            return false; 
        }

        FILE* fp = fopen("cut_plc_faces.off", "w");
        fprintf(fp, "OFF\n%zu %zu 0\n", vertices.size(), faces.size());
        for(uint32_t i=0; i<vertices.size(); i++) {
            assert(vertices[i]->is3D()); // DEBUG
            double x, y, z;
            vertices[i]->getApproxXYZCoordinates(x, y, z);
            fprintf(fp, "%f %f %f\n", x, y, z);
        }
        
        for(const CHAMface& f : faces) {
            fprintf(fp, "%zu ", f.bounding_edges.size());
            std::vector<uint32_t> fv;
            get_face_vertices(f, fv);
            for(size_t i=0; i<fv.size(); i++) fprintf(fp, "%u ", fv[i]);
            fprintf(fp, "\n");
        }
        fclose(fp);
        return true;
    }

    // computes a trinagulation of the chamfered PLC (without modifing it) 
    // and saves into the result in to an off file
    bool saveTriFaces() const {

        for(const CHAMface& f : faces)if( !f.is_simply_connected ){  
            std::cout<<"PLCc - not simply connected face founded (off file not generated)\n";
            return false; 
        }

        std::vector<uint32_t> out_tri;
        get_triangles(out_tri);
        assert( (out_tri.size() % 3) == 0 );

        uint32_t nfaces = (uint32_t)out_tri.size() / 3;

        FILE* fp = fopen("cut_plc_faces_triangles.off", "w");
        fprintf(fp, "OFF\n%zu %u 0\n", vertices.size(), nfaces);
        for(uint32_t i=0; i<vertices.size(); i++) {
            assert(vertices[i]->is3D()); // DEBUG
            double x, y, z;
            vertices[i]->getApproxXYZCoordinates(x, y, z);
            fprintf(fp, "%f %f %f\n", x, y, z);
        }
        
        for(size_t i=0; i<out_tri.size()/3; i++) {
            std::vector<uint32_t> v;
            v.assign(out_tri.begin() + 3*i, out_tri.begin() +3*i +3);
            fprintf(fp, "3 %u %u %u \n", v[0], v[1], v[2]);
        }
        fclose(fp);
        return true;
    }

    // save the triangles whose vertices indices are triple of consecutive integers
    // of save_tri vector into a an .off file
    bool saveTriangles(const std::vector<uint32_t>& save_tri, const char* title) const {

        assert( (save_tri.size() % 3) == 0 );

        std::vector<uint32_t> used_vertex(vertices.size(), UINT32_MAX);
	    for (size_t i = 0; i < save_tri.size(); i++)    used_vertex[save_tri[i]] = 1;

	    uint32_t idx = 0;
	    for (size_t i = 0; i < used_vertex.size(); i++) if (used_vertex[i] != UINT32_MAX) {
		    used_vertex[i] = idx++;
	    }

        uint32_t nfaces = (uint32_t)save_tri.size() / 3;
        FILE* fp = fopen(title, "w");
        fprintf(fp, "OFF\n%u %u 0\n", idx, nfaces);

        for(uint32_t i=0; i<vertices.size(); i++) if(used_vertex[i] != UINT32_MAX){
            assert(vertices[i]->is3D()); // DEBUG
            double x, y, z;
            vertices[i]->getApproxXYZCoordinates(x, y, z);
            fprintf(fp, "%f %f %f\n", x, y, z);
        }

        for(size_t i=0; i<save_tri.size()/3; i++) {
            std::vector<uint32_t> v;
            v.assign(save_tri.begin() + 3*i, save_tri.begin() +3*i +3);
            fprintf(fp, "3 %u %u %u \n", used_vertex[v[0]], used_vertex[v[1]], used_vertex[v[2]]);
        }
        
        fclose(fp);
        return true;
    }

    // do not work for non-simpy connected faces
    bool saveFace(const CHAMface& f, const char* title) const {

        if( !f.is_simply_connected ){  
            std::cout<<"PLCc - saveFace() - not simply connected face (off file not generated)\n";
            return false; 
        }

        std::vector<uint32_t> fv;
        get_face_vertices(f, fv);

        FILE* fp = fopen(title, "w");
        fprintf(fp, "OFF\n%zu 1 0\n", fv.size());
        for(uint32_t i=0; i<fv.size(); i++) {
            double x, y, z;
            vertices[fv[i]]->getApproxXYZCoordinates(x, y, z);
            fprintf(fp, "%f %f %f\n", x, y, z);
        }
        
        fprintf(fp, "%zu ", f.bounding_edges.size());

        for(uint32_t i=0; i<fv.size(); i++) fprintf(fp, "%u ", i);
        fprintf(fp, "\n");
        fclose(fp);

        return true;
    }

    void saveFaceVertices(const std::vector<uint32_t> fv, const char* title) const {

        FILE* fp = fopen(title, "w");
        fprintf(fp, "OFF\n%zu 1 0\n", fv.size());
        for(uint32_t i=0; i<fv.size(); i++) {
            double x, y, z;
            vertices[fv[i]]->getApproxXYZCoordinates(x, y, z);
            fprintf(fp, "%f %f %f\n", x, y, z);
        }
        
        fprintf(fp, "%zu ", fv.size());

        for(uint32_t i=0; i<fv.size(); i++) fprintf(fp, "%u ", i);
        fprintf(fp, "\n");
        fclose(fp);
    }

    void saveTriangle(uint32_t v1, uint32_t v2, uint32_t v3, const char* title) const {
        FILE* fp = fopen(title, "w");
        fprintf(fp, "OFF\n3 1 0\n");
        double x, y, z;
        vertices[v1]->getApproxXYZCoordinates(x, y, z);
        fprintf(fp, "%f %f %f\n", x, y, z);
        vertices[v2]->getApproxXYZCoordinates(x, y, z);
        fprintf(fp, "%f %f %f\n", x, y, z);
        vertices[v3]->getApproxXYZCoordinates(x, y, z);
        fprintf(fp, "%f %f %f\n", x, y, z);
        fprintf(fp, "3 0 1 2\n");
        fclose(fp);
    }

    void saveTriangle(const double* p1c, const double* p2c, const double* p3c, const char* title) const {
        FILE* fp = fopen(title, "w");
        fprintf(fp, "OFF\n3 1 0\n");
        fprintf(fp, "%f %f %f\n", p1c[0], p1c[1], p1c[2]);
        fprintf(fp, "%f %f %f\n", p2c[0], p2c[1], p2c[2]);
        fprintf(fp, "%f %f %f\n", p3c[0], p3c[1], p3c[2]);
        fprintf(fp, "3 0 1 2\n");
        fclose(fp);
    }

    // do not work if non-simpy connected faces are present
    void saveEdgeIncFaces(const CHAMedge& e) const {

        for(uint32_t fi : e.inc_face)if( !faces[fi].is_simply_connected ){  
            std::cout<<"PLCc - saveEdgeIncFaces() - not simply connected face founded (off file not generated)\n";
            return; 
        }

        FILE* fp = fopen("edge_inc_faces.off", "w");
        fprintf(fp, "OFF\n%zu %zu 0\n", vertices.size(), e.inc_face.size());
        for(uint32_t i=0; i<vertices.size(); i++) {
            double x, y, z;
            vertices[i]->getApproxXYZCoordinates(x, y, z);
            fprintf(fp, "%f %f %f\n", x, y, z);
        }
        
        for (uint32_t fi : e.inc_face) {
            const CHAMface& f = faces[ fi ];
            fprintf(fp, "%zu ", f.bounding_edges.size());

            std::vector<uint32_t> fv;
            get_face_vertices(f, fv);
            for(size_t i=0; i<fv.size(); i++) fprintf(fp, "%u ", fv[i]);
            
            fprintf(fp, "\n");
        }
        fclose(fp);
    }

    void save_rebuilded_input_after_chamfering(const std::vector<uint32_t>& tri_cham_faces, const char* title){
        std::vector<uint32_t> compl_tri;
        get_complementar_tri(tri_cham_faces, compl_tri);
		compl_tri.insert(compl_tri.end(), tri_cham_faces.begin(), tri_cham_faces.end());
		saveTriangles(compl_tri, title);
    }


    // Debug

    inline void print_vertex(uint32_t vi) const { std::cout<<"vertex["<<vi<<"] = "<< *vertices[vi] <<"\n"; }
    inline void print_vertex_info(uint32_t vi) const { 
        const pointType* v = vertices[vi];
        std::cout<<"vertex["<<vi<<"] = "<< *v
                 <<" - type ";
        if(v->isExplicit3D()) std::cout<<"EXPLIC3D ";
        else if(v->isLNC()){ std::cout<<"LNC (P = "<< v->toLNC().P() <<", Q = "<< v->toLNC().Q() <<", t = "<<v->toLNC().T()<<" )"; }
        else if(v->isBPT()){ std::cout<<"BPT (P = "<< v->toBPT().P() <<", Q = "<< v->toBPT().Q() <<", R = "<< v->toBPT().R() <<", u = "<<v->toBPT().U()<<", v = "<<v->toBPT().V()<<" )"; }
        else std::cout<<"UNKNOWN ";
        std::cout<<"\n"; 
    }
    inline void print_edge(size_t ei) const { std::cout<<"edge["<<(uint32_t)ei<<"] = <"<<edges[ei].ep[0]<<","<<edges[ei].ep[1]<<">\n"; }
    inline void print_input_vt(uint32_t vi) const { 
        std::cout<<"input vt["<<vi<<"]:\n";
        for(size_t i=0; i<input_vt[vi].size(); i++){
            const uint32_t* ti = plc.triangle_vertices.data() + input_vt[vi][i]*3;
            std::cout<<*ti<<" "<<*(ti+1)<<" "<<*(ti+2)<<std::endl;
        }
        std::cout<<std::endl;
     }
    void print_edge_and_inc_face(uint32_t ei) const { 
        std::cout<<"edge["<<ei<<"] = <"<<edges[ei].ep[0]<<","<<edges[ei].ep[1]<<"> "
                 <<" inc_face: ";
        for(uint32_t fi : edges[ei].inc_face) std::cout<<fi<<" ";
        std::cout<<"\n"; 
    }
    void print_face(uint32_t fi) const {
        std::cout<<"face["<<fi<<"] : ("<< faces[fi].bounding_edges.size() <<" bound-edges) \n";
        for(size_t i = 0; i < faces[fi].bounding_edges.size(); i++){  
            std::cout<<"edge["<< faces[fi].bounding_edges[i] <<"]\n"; 
        }
        std::cout<<"\n";
    }
    void print_face_edges(uint32_t fi) const {
        std::cout<<"face["<<fi<<"] : ("<< faces[fi].bounding_edges.size() <<" bound-edges) \n";
        for(uint32_t ei : faces[fi].bounding_edges)
            if(ei!=UINT32_MAX){  print_edge(ei); }
            else std::cout<<"JUNK\n";
        std::cout<<"\n";
    }

    bool vertex_have_sense(uint32_t v) const {
        bool error = false;
        if( v >= vertices.size()){
            error = true;
            std::cout<<"ERROR PLCc - vertex_have_sense(): "
                     <<"vertex index "<<v<<" greater than vertices size ("<< vertices.size() <<")\n";
        }
        if( !vertices[v]->is3D() ){
            error = true;
            std::cout<<"ERROR PLCc - vertex_have_sense(): "
                     <<"vertex["<<v<<"] is not 3D (inpute vertices are first "<< n_in_vrts <<")\n";

        }
        return !error;
    }

    // during vertex chamfering
    void check_bridge(uint32_t v1, uint32_t v2, uint32_t v3, uint32_t v4) const{
        if(!vertices[v1]->isLNC() || !vertices[v4]->isLNC() ||
           !vertices[v2]->isBPT() || !vertices[v3]->isBPT() ){
                uint32_t v[] = {v1,v2,v3,v4};
                for(size_t i = 0; i < 4; i++){
                    std::cout<<"v"<<i<<" -> ";
                    if(vertices[v[i]]->isLNC()) std::cout<<"LNC\n";
                    else if(vertices[v[i]]->isBPT()) std::cout<<"BPT\n";
                    else if(vertices[v[i]]->isExplicit3D()) std::cout<<"EXP3D\n";
                    else std::cout<<"OTHER\n";
                }
                // ip_error("[check_bridge] ERROR invalid bridge vertices\n");
                std::cout<<"[cham.h - check_bridge()] ERROR invalid bridge vertices\n"; exit(1);
           }
    }

    // during simplification
    void print_all_bridge_edges() const {
        std::cout<<"\nBRIDGE EDGES:\n";
        for(uint32_t ei = 0; ei < edges.size(); ei++) if(edges[ei].loc_face_bridge_id != NO_BRIDGE) { 
            const CHAMedge& e = edges[ei];
            std::cout<<"\n"; 
            std::cout<<"<";
            if( vertices[e.ep[0]]->isLNC() ) std::cout<<"LNC, ";
            else if( vertices[e.ep[0]]->isBPT() ) std::cout<<"BPT, ";
            else if( vertices[e.ep[0]]->isExplicit3D() ) std::cout<<"!EXPL3D!, ";
            else std::cout<<"!OTHER!, ";
            if( vertices[e.ep[1]]->isLNC() ) std::cout<<"LNC> -> ";
            else if( vertices[e.ep[1]]->isBPT() ) std::cout<<"BPT> -> ";
            else if( vertices[e.ep[1]]->isExplicit3D() ) std::cout<<"!EXPL3D!> -> ";
            else std::cout<<"!OTHER!> -> ";
            print_edge(ei);
            std::cout<<"lfb_id = "; if(e.loc_face_bridge_id!=NO_BRIDGE) std::cout<<e.loc_face_bridge_id<<"\n"; else std::cout<<"NO_BRIGE\n";
            std::cout<<"n_inc_faces = "<<e.inc_face.size()<<"\n";
        }
        std::cout<<std::endl;

    }

    void check_bridge_edges() const {
        for(const CHAMedge& e : edges) if(e.loc_face_bridge_id != NO_BRIDGE) {
            assert( e.inc_face.size()==1);
            assert( vertices[ e.ep[0] ]->isLNC() || vertices[ e.ep[0] ]->isBPT() );
            assert( vertices[ e.ep[1] ]->isLNC() || vertices[ e.ep[1] ]->isBPT() );
            assert( !( vertices[ e.ep[0] ]->isLNC() && vertices[ e.ep[1] ]->isLNC()) );
        }
    }
    
    // check that each acute edge has both acute endpoints: use it only at the end of initialize()
    bool acute_edges_have_acute_ep() const {
        bool error = false;
        for(uint32_t ei=0; ei<edges.size(); ei++) if(edges[ei].isAcute()){ 
            const CHAMedge& e = edges[ei];
            if( !is_acute_vrt(e.ep[0]) || !is_acute_vrt(e.ep[1]) ){
                error = true;
                uint32_t n = (uint32_t) (!is_acute_vrt(e.ep[0])) + (uint32_t) (!is_acute_vrt(e.ep[1]));
                std::cout<<"ERROR PLCc - check_acute_coherence(): "
                         <<"acute edge["<<ei<<"] has "<<n<<" not-acute endpoints\n";
                
                //for(uint32_t fi : e.inc_face) print_face_edges( fi );
                //saveEdgeIncFaces(e); assert(is_acute_vrt(e.ep[0]) && is_acute_vrt(e.ep[1]));
            }
        }
        return(!error);
    }

    bool edge_is_coherent(uint32_t ei) const {
        const CHAMedge& e = edges[ei];
        bool error = false;

        if( ei >= edges.size()){
            error = true;
            std::cout<<"ERROR PLCc - edge_is_coherent(): "
                     <<"edge index "<<ei<<" greater than edges size ("<< edges.size() <<")\n";
        }
        
        // check endpoints
        if( !vertex_have_sense(e.ep[0]) || !vertex_have_sense(e.ep[1]) ){ error = true; }

        #ifndef FP_CHAMFERING
        
        if(e.oep[0]!=e.ep[0] || e.oep[1]!=e.ep[1]){

            if(e.oep[0]!=e.ep[0] && !vertices[e.ep[0]]->isBPT()){
                if(!vCollinear(e.oep[0], e.ep[0], e.ep[1])){
                    std::cout<<"ERROR PLCc - edge_is_coherent(): "
                             <<"edge["<<ei<<"] oep0, ep0, ep1 NOT collinear\n";
                    error = true;
                }

                if(!vPointInInnerSegment(e.ep[0], e.oep[0], e.ep[1])){
                    std::cout<<"ERROR PLCc - edge_is_coherent(): "
                             <<"edge["<<ei<<"] ep0 outside <oep0,ep1>\n";
                    error = true;
                }
            }

            if(e.oep[1]!=e.ep[1] && !vertices[e.ep[1]]->isBPT()){ 
                if(!vCollinear(e.ep[0], e.ep[1], e.oep[1])){
                    std::cout<<"ERROR PLCc - edge_is_coherent(): "
                             <<"edge["<<ei<<"] ep0, ep1, oep1 NOT collinear\n";
                    error = true;
                }

                if(!vPointInInnerSegment(e.ep[1], e.ep[0], e.oep[1])){
                    std::cout<<"ERROR PLCc - edge_is_coherent(): "
                             <<"edge["<<ei<<"] ep1 outside <ep0,oep1>\n";
                    error = true;
                }

            }
        }
        #endif

        return !error;
    }

    bool face_boundary_overlaps(uint32_t fi) const {
        const CHAMface& f = faces[fi];
        const std::vector<uint32_t>& fbe = f.bounding_edges;
        for(size_t i = 0; i < fbe.size()-1; i++){
            const CHAMedge& ei = edges[ fbe[i] ];
            for(size_t j = i+1; j < fbe.size(); j++){
                const CHAMedge& ej = edges[ fbe[j] ];

                if( vInnerSegmentCross(ei.ep[0], ei.ep[1], ej.ep[0], ej.ep[1]) ){
                    std::cout<<"ERROR PLCc - face_boundary_overlaps(): "
                        <<"edge["<< fbe[i] <<"] and edge["<< fbe[j] <<" intersects\n";
                    print_face_edges(fi);
                    saveFace(f, "face.off");
                    return false;
                }

            }
                        
        }

        return true;
    }

    bool face_boundary_have_sense(uint32_t fi) const {
        bool error = false;
        const CHAMface& f = faces[fi];
        const std::vector<uint32_t>& fbe = f.bounding_edges;

        // check tath each edge has at two neighbour (with which share an endpoint)
        for(size_t it=0; it < fbe.size(); ++it){
            size_t n_it = it;  f.advance_on_bnd(n_it);
            size_t p_it = it;  f.reverse_on_bnd(p_it);

            const CHAMedge& e = edges[ fbe[it] ];
            const CHAMedge& ne = edges[ fbe[n_it] ];
            const CHAMedge& pe = edges[ fbe[p_it] ];

            if(!e.has_commonVertex(ne) || !e.has_commonVertex(pe) ||
                e.coincident(ne) || e.coincident(pe)           ){
                    std::cout<<"ERROR PLCc - face_boundary_have_sense(): "
                             <<"edge["<<fbe[it]<<"] have problems with neighbours\n";
                    print_face_edges(fi);
                    error = true;

            }
        }

        // check the correct orientation wrt neighbours on boundary
        #ifndef FP_CHAMFERINF
        if( !face_boundary_overlaps(fi) ) error = true;
        #endif

        // check that each enpoint compares exactly two times on the bounday
        std::vector<uint32_t> loc_mark_vrts(vertices.size(), 0);
        for(uint32_t ei : f.bounding_edges){
            loc_mark_vrts[edges[ei].ep[0]]++;
            loc_mark_vrts[edges[ei].ep[1]]++;
        }

        for(uint32_t ei : f.bounding_edges) for(uint32_t k=0; k<2; k++) {
            uint32_t vi = edges[ei].ep[k];
            if(loc_mark_vrts[vi]!=2){
                std::cout<<"ERROR PLCc - face_boundary_have_sense(): "
                    <<"endpoint "<<vi<<" compares "<<loc_mark_vrts[vi]<<" times\n";
                error = true;
            }
        }

        
            
        return !error;
    }

    bool EF_relation_is_choerent(uint32_t ei) const{
        const CHAMedge& e = edges[ei];
        for(uint32_t fi : e.inc_face) if(!faces[fi].hasBndEdge(ei)){ 
            std::cout<<"ERROR PLCc - check_EF_relation(): "
                     <<"edge["<<ei<<"] has inc_face face["<<fi<<"] "
                     <<"but do not appear on its boundary\n";
            return false; 
        }
        return true;
    }

    bool FE_relation_is_choerent(uint32_t fi) const{
        const CHAMface& f = faces[fi];
        for(uint32_t ei : f.bounding_edges){ 
            const CHAMedge& e = edges[ei];
            if(std::find(e.inc_face.begin(), e.inc_face.end(), fi) == e.inc_face.end()){ 
                std::cout<<"ERROR PLCc - check_FE_relation(): "
                        <<"face["<<fi<<"] has bnd_edge edge["<<ei<<"] "
                        <<"but do not appear as its incident face\n";
                return false; 
            }
        }
        return true;
    }

    bool checkup() const {
        bool error = false;

        // check edges
        #ifndef FP_CHAMFERING
        for(uint32_t ei=0; ei<edges.size(); ei++) if(!edge_is_coherent(ei)) error = true;
        #endif
        
        // check faces
        for(uint32_t fi=0; fi<faces.size(); fi++) if(!face_boundary_have_sense(fi)) error = true;

        // check connectivity
        for(uint32_t ei=0; ei<edges.size(); ei++) if(!EF_relation_is_choerent(ei)) error = true;
        for(uint32_t fi=0; fi<faces.size(); fi++) if(!FE_relation_is_choerent(fi)) error = true;

        return !error;
    }

    bool analyze_acute_angle(const pointType* pl, const pointType* p, const pointType* pr, uint32_t v, uint32_t fi) const {

        const double cos_toll = 0.0157; // approximation of cos(pi/2 * 99/100) ~= cos(89,1°)

        vector3d Opl(pl), Opr(pr), Op(p);

        const double comp_dp = abs((Opl-Op).dot(Opr-Op));
        const double sq_ll = Opl.dist_sq(Op); 
        const double sq_lr = Opr.dist_sq(Op); 
        const double norm_prod = sqrt( sq_ll * sq_lr );

        if( comp_dp < cos_toll*norm_prod ){
            #ifdef PLCC_VERBOSE_DEBUG
            std::cout<<" vertex "<< v <<" (type "<<vertices[v]->getType()<<") on face["<<fi<<"]:\n";
            std::cout<<"[ analyze_acute_angle() ] WARNING: cos(alpha) = "<< comp_dp / norm_prod <<"\n";
            #endif
            return false;
        }

        std::cout<<" \nvertex "<< v <<" (type "<<vertices[v]->getType()<<") on face["<<fi<<"]:\n";
        std::cout<<"[ analyze_acute_angle() ] ERROR: cos(alpha) = "<< comp_dp / norm_prod <<"\n";
        std::cout<<"[ analyze_acute_angle() ] comp_dotprod = "<< comp_dp <<"\n";
        std::cout<<"[ analyze_acute_angle() ] norm_prod = "<< norm_prod <<"\n";
        std::cout<<"[ analyze_acute_angle() ] cent node (v) coords = "<< Op <<"\n";
        std::cout<<"[ analyze_acute_angle() ] vertices incident at v: "<<Opl<<" , "<<Opr<<"\n";
        uint32_t intri = faces[fi].triangle;
        const uint32_t* intriv = plc.triangle_vertices.data() + intri*3;
        std::cout<<"[ analyze_acute_angle() ] original triangle "<<intri<<": <"<< *intriv <<","<< *(intriv+1) << ","<< *(intriv+2) <<">\n";
        Opl = (Op+Opl)*0.5;
        Opr = (Op+Opr)*0.5;
        saveTriangle(Opl.c,Op.c,Opr.c,"acute_angle.off");
        return true;

    }

    // use it only at the end of chamfering algorithm
    bool check_acuteness() const {

        // check that each couple of faces incident at an edges do not form 
        // an acute dihedral angle
        // DO NOT WORK CAUSE OF NOT-COVEX FACES 
        // for(uint32_t ei=0; ei<edges.size(); ei++){

        //     const CHAMedge& e = edges[ei];

        //     if(e.inc_face.size() == 1) continue;
        //     if(e.inc_face.size() > 4){ 
        //         std::cout<<"[ check_acuteness() ] ERROR: edges["<<ei<<"] has "
        //                  << e.inc_face.size() <<" incident faces, thus is acute.\n";
        //         return false;
        //     }

        //     uint32_t e0 = e.ep[0], e1 = e.ep[1];
        //     uint32_t u, v;
        //     const pointType* e0_pt = vertices[e0];
        //     const pointType* e1_pt = vertices[e1];

        //     for(size_t i=1; i<e.inc_face.size(); i++){
        //         for(size_t j=0; j<i; j++){
        //             u = get_not_collinear_vertex_on_face(e0, e1, e.inc_face[i]);
        //             v = get_not_collinear_vertex_on_face(e0, e1, e.inc_face[j]);

        //             const vector3d Oe0( e0_pt );
        //             const vector3d e0e1( vector3d(e1_pt) - Oe0 );
        //             const vector3d nu = e0e1.cross( vector3d( vertices[u] ) - Oe0 );
        //             const vector3d nv = e0e1.cross( vector3d( vertices[v] ) - Oe0 );
        //             const vector3d Onu( Oe0 + nu );
        //             const vector3d Onv( Oe0 + nv );

        //             const explicitPoint3D p( Onu.c[0], Onu.c[1], Onu.c[2] ) ;
        //             const explicitPoint3D q( Onv.c[0], Onv.c[1], Onv.c[2] ) ;

        //             if( pointType::dotProductSign3D(p, q, *e0_pt) > 0){
        //                 std::cout<<"\ns_dp = "<<pointType::dotProductSign3D(p, q, *e0_pt)
        //                          <<" for face["<<e.inc_face[i]<<"] and face["<<e.inc_face[j]<<"] incident at ";
        //                 print_edge(ei);
        //                 if( analyze_acute_angle(&p, e0_pt, &q) ){ 
        //                     const vector3d Oe1( e1_pt );
        //                     const double pc[3] = {p.X(), p.Y(), p.Z()};
        //                     const double qc[3] = {q.X(), q.Y(), q.Z()};
        //                     saveTriangle(pc, Oe0.c, Oe1.c, "norm1.off");
        //                     saveTriangle(qc, Oe0.c, Oe1.c, "norm2.off");
        //                     if( !saveFace(faces[ e.inc_face[i] ], "face1.off") )
        //                         saveTriangle(e0,e1,u,"tri1.off");
        //                     if( !saveFace(faces[ e.inc_face[j] ], "face2.off"))
        //                         saveTriangle(e0,e1,v,"tri2.off");
        //                     return false;

        //                 }
        //             }
        //         }
        //     }

        // }

        // check that each face as polygon has not acute angles
        for(uint32_t fi=0; fi<faces.size(); fi++){
            const CHAMface& f = faces[fi];
            const std::vector<uint32_t>& fbnd = f.bounding_edges;
            size_t it = 0;
            while(it < fbnd.size()){
                size_t nit = it; f.advance_on_bnd(nit);
                uint32_t v, vl, vr;
                const CHAMedge& e1 = edges[ f.bounding_edges[it] ];
                const CHAMedge& e2 = edges[ f.bounding_edges[nit] ];
                v = e1.commonVertex(e2);
                vl = e1.notCommonVertex(e2);
                vr = e2.notCommonVertex(e1);
                //if( pointType::dotProductSign3D( *vertices[vl], *vertices[vr], *vertices[v]) > 0 ){ 
                if( !findIF_flat_edge(e1) && 
                    !findIF_flat_edge(e2) &&
                    isAcuteAngle( vertices[vl], vertices[v], vertices[vr]) ){ 
                    if( analyze_acute_angle(vertices[vl], vertices[v], vertices[vr],v,fi) ){ 
                        saveFace(f,"acute_face.off");
                        return false;
                    }
                }
                ++it;
            }
        }

        // check that there are no acute vertices (considering all incident elements)


        return true;
    }

    #ifdef PLCC_DEBUG
    // Debug mode functions

    bool duplicated_vertices() const{
        bool duplicated = false;
        std::vector< std::pair<vector3d,uint32_t> > vrts_cpy;
        for(size_t i=0; i<vertices.size(); i++){ 
            vrts_cpy.push_back( std::pair<vector3d,uint32_t> (vector3d(vertices[i]),i) ); 
        }
        std::sort(vrts_cpy.begin(), vrts_cpy.end(),
                  [](const std::pair<vector3d,uint32_t> &a, const std::pair<vector3d,uint32_t> &b){ return a.first < b.first;} );
        for(size_t i=1; i<vrts_cpy.size(); i++){ 
            if(vrts_cpy[i-1].first.c[0] == vrts_cpy[i].first.c[0] &&
            vrts_cpy[i-1].first.c[1] == vrts_cpy[i].first.c[1] &&
            vrts_cpy[i-1].first.c[2] == vrts_cpy[i].first.c[2]     ){
                uint32_t v1 = vrts_cpy[i-1].second;
                uint32_t v2 = vrts_cpy[i].second;
                std::cout<<"ERROR: duplicated vertex["<< v1 <<"] = vertex["<< v2 <<"]\n"; 
                print_vertex_info(v1);
                print_vertex_info(v2);
                const pointType *vv1 = vertices[v1];
                const pointType *vv2 = vertices[v2];
                if(vv1->isLNC() && vv2->isLNC()){  
                    if(vv1->toLNC().P() == vv2->toLNC().P() && 
                       vv1->toLNC().Q() == vv2->toLNC().Q() && 
                       vv1->toLNC().T() == vv2->toLNC().T()   )
                        std::cout<<"Exact Versions Are Equal\n";
                }
                if(vv1->isBPT() && vv2->isBPT()){  
                    if(vv1->toBPT().P() == vv2->toBPT().P() && 
                       vv1->toBPT().Q() == vv2->toBPT().Q() && 
                       vv1->toBPT().R() == vv2->toBPT().R() && 
                       vv1->toBPT().U() == vv2->toBPT().U() &&
                       vv1->toBPT().V() == vv2->toBPT().V()   )
                        std::cout<<"Exact Versions Are Equal\n";
                }
                
                duplicated = true;
            }
        }
        return duplicated;
    }

    bool duplicated_edges() const{
        bool duplicated = false;
        std::vector< std::vector<uint32_t> > edge_cpy;
        for(uint32_t i=0; i<edges.size(); i++){ 
            uint32_t ep0=edges[i].ep[0], ep1=edges[i].ep[1];
            if(ep0>ep1) std::swap(ep0,ep1);
            edge_cpy.push_back( std::vector<uint32_t> ({ep0, ep1, i}) ); 
        }
        std::sort(edge_cpy.begin(), edge_cpy.end(),
                  [](const std::vector<uint32_t> &a, const std::vector<uint32_t> &b){ 
                                        return (a[0] < b[0] || (a[0]==b[0] && a[1]<b[1]));
                                        } 
                  );
        for(size_t i=1; i<edge_cpy.size(); i++){ 
            if(edge_cpy[i-1][0] == edge_cpy[i][0] &&
               edge_cpy[i-1][1] == edge_cpy[i][1]    ){
                std::cout<<"ERROR: duplicated edge["<< edge_cpy[i-1][2] <<"] = edge["<< edge_cpy[i][2] <<"]\n"; 
                duplicated = true;
            }
        }
        return duplicated;
    }

    void checkEdges_beforeJunkDeletion() const {
        for(const CHAMedge& e : edges)if(!e.isIsolated()){
            if(e.isJunk()) assert( is_acute_vrt(e.ep[0]) || is_acute_vrt(e.ep[1]) );
            else assert( !is_acute_vrt(e.ep[0]) && !is_acute_vrt(e.ep[1]) );
        }
    }

    void checkFacesBND_afterJunkDeletion() const {
        for(uint32_t fi=0; fi<faces.size(); fi++){
            const CHAMface& f = faces[fi];
            for(uint32_t ei : f.bounding_edges){
                const CHAMedge& e = edges[ei];
                assert(!is_acute_vrt(e.ep[0]));
                assert(!is_acute_vrt(e.ep[1]));
            }
        }
    }

    //
    void report_vons_fail(uint32_t v0, uint32_t v1){
        string type_ep0, type_ep1;

        if(vertices[v0]->isExplicit3D()) type_ep0 = "explicit";
        else if(vertices[v0]->isLNC()) type_ep0 = "LNC";
        else if(vertices[v0]->isBPT()) type_ep0 = "BTP";
        else type_ep0 = "unknown";

        if(vertices[v1]->isExplicit3D()) type_ep1 = "explicit";
        else if(vertices[v1]->isLNC()) type_ep1 = "LNC";
        else if(vertices[v1]->isBPT()) type_ep1 = "BTP";
        else type_ep1 = "unknown";

        std::cout<<"[new_vrt_on_segment()] ERROR: "
                 <<"ep0 = "<<v0<<" is " << type_ep0 << " and "
                 <<"ep1 = "<<v1<<" is " << type_ep1 <<": this shold never happen.\n";
    }

    #endif



    #ifdef PLCC_VERBOSE_DEBUG
    // Verbose mode functions

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

    #endif

    #ifdef PLCC_VERBOSE_DEBUG_LEV1

    // More-verbose mode functions
    inline void disp_allFlatEdges() const { std::cout<<"FLAT EDGES:\n"; for(uint32_t i=0; i<edges.size(); i++) if(edges[i].isFlat()) print_edge(i); std::cout<<endl;}
    inline void disp_allAcuteEdges() const { std::cout<<"ACUTE EDGES:\n"; for(uint32_t i=0; i<edges.size(); i++) if(edges[i].isAcute()) print_edge(i); std::cout<<endl;}
    inline void disp_allAcuteVertices() const { std::cout<<"ACUTE VRTS:\n"; for(uint32_t i=0; i<vertices.size(); i++) if(is_acute_vrt(i)) print_vertex(i); std::cout<<endl;}
    #else
    //inline void disp_allFlatEdges() const { std::cout<<"FLAT EDGES:\n"; for(uint32_t i=0; i<edges.size(); i++) if(edges[i].isFlat()) print_edge(i); std::cout<<endl;}

    #endif


};

#include "cham.cpp"

#endif // _CHAM_
