#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <set>

// ---------- //
// CHAM FACES //
// ---------- //

void CHAMface::make_first(uint32_t ei){
    u32vect_iter it = std::find(bounding_edges.begin(), bounding_edges.end(), ei);
    if(it!=bounding_edges.end()) make_first(it);
}

void CHAMface::make_last(u32vect_iter& it){
    // To make the edge e (pointed by it) the last of bounding_edges 
    // we rotate bounding_edges in such a way that e+1 (the succesive on boundary of e) 
    // occupy the first position.
    advance_on_bnd(it);
    std::rotate( bounding_edges.begin(), it, bounding_edges.end() );
}

void CHAMface::make_last(uint32_t ei){
    u32vect_iter it = std::find(bounding_edges.begin(), bounding_edges.end(), ei);
    if(it!=bounding_edges.end()) make_last(it);
}

void CHAMface::replaceEdge_11(uint32_t old_e, uint32_t new_e){
    // NOTE: replace all istances of old_e with new_e
    //       but we assume that old_e is unique on bounding edges.
    std::replace(bounding_edges.begin(), bounding_edges.end(), old_e, new_e);
    //for(uint32_t& e : bounding_edges) if(e==old_e){ e=new_e; break; }
}

// ------ //
// ANGLES //
// ------ //

// Let be tau the plane for the triangle <q,r,s>.
// Returns TRUE if 
// segment qp forms an acute angle at q with the triangle <q,r,s> AND
// the projection of qp on tau is inside the sector limited by half-straight-lines qr and qs
// containing the triangle <q,r,s>
bool isAcuteAngle(const pointType* p, const pointType* q, const pointType* r, const pointType* s) {
    // The following check is based on the variatinal approach described in
    // https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
    
    vector3d Op(p), Oq(q), Or(r), Os(s), w(Oq-Op);
    double a = Or.dot(Or), b = Or.dot(Os), c = Os.dot(Os);
    double d = Or.dot(w), e = Os.dot(w);
    double u = b * e - c * d;
    double v = b * d - a * e;
    return (u > 0 && v > 0);
}

// ------------------ //
// IMPLICIT RELATIONS //
// ------------------ //

// Supports only simply-connected faces (each vertex compares at most 1 time on face boundary).
// Fills fv with the indices of face f vertices.
void PLCc::get_face_vertices(const CHAMface& f, std::vector<uint32_t>& fv) const {

    const std::vector<uint32_t>& fbe = f.bounding_edges;

    fv.resize( fbe.size() );
    uint32_t v = UINT32_MAX;
    for(size_t i = 0; i<fbe.size()-1; i++){
        const CHAMedge& e1 = edges[ fbe[i] ];
        const CHAMedge& e2 = edges[ fbe[i+1] ];

        #ifdef PLCC_DEBUG
        assert( e1.has_commonVertex( e2 ) );
        #endif
        
        fv[i] = e1.commonVertex( e2 );
    }

    #ifdef PLCC_DEBUG
    assert( edges[ fbe.back() ].has_commonVertex( edges[ fbe.front() ] ) );
    #endif
        
    fv.back() = edges[ fbe.back() ].commonVertex( edges[ fbe.front() ] );
}

// ------------------------- //
// OREDERING OF PLCC VECTORS //
// ------------------------- //

void PLCc::swap_edges(uint32_t e1, uint32_t e2){

    // update boundary edges of faces incident at edges[e1]
    for(size_t i=0; i<edges[e1].inc_face.size(); i++)
        faces[ edges[e1].inc_face[i] ].replaceEdge_11(e1, e2);

    // update boundary edges of faces incident at edges[e2]
    for(size_t i=0; i<edges[e2].inc_face.size(); i++)
        faces[edges[e2].inc_face[i]].replaceEdge_11(e2, e1);

    std::swap(edges[e1], edges[e2]); // swap edges elements
    std::swap(mark_edges[e1], mark_edges[e2]);
}

uint32_t PLCc::move_back_isolated_edges(){

    uint32_t last = ( (uint32_t) edges.size() ) -1;
    for(uint32_t i=0; i<edges.size(); i++) if(edges[i].isIsolated()){ 
        while( edges[last].isIsolated() ) last--;
        if(last > i){ swap_edges(i, last); }
    }
    while( edges[last].isIsolated() ) last--;
    return ++last;

}

// ------------------- //
// PLCC INITIALIZATION //
// ------------------- //

// returns TRUE if triangle u = <u[0],u[1],u[2]> is flipped
// wrt triangle v = <v[0],v[1],v[2]>
static inline bool are_flipped_triangles(uint32_t u[], uint32_t v[]){
    return( (v[0]==u[1] && v[1]==u[0] && v[2]==u[2]) ||
            (v[0]==u[2] && v[1]==u[1] && v[2]==u[0]) ||
            (v[0]==u[0] && v[1]==u[2] && v[2]==u[1])    );
}

void PLCc::orient_initial_triface_bnd(CHAMface& f){
    const uint32_t* t = plc.triangle_vertices.data() + (f.triangle*3);
    uint32_t v[] = {*t, *(t+1), *(t+2)};
    const CHAMedge& e0 = edges[ f.bounding_edges[0] ];
    const CHAMedge& e1 = edges[ f.bounding_edges[1] ];
    const CHAMedge& e2 = edges[ f.bounding_edges[2] ];
    uint32_t u[] = {e0.commonVertex(e1), e1.commonVertex(e2), e2.commonVertex(e0)};
    if( are_flipped_triangles(u,v) ) std::swap(f.bounding_edges[0], f.bounding_edges[2]);
}

// Removes duplicated pre-edges (treansorming half-edges into edges) 
void PLCc::mergePreEdges(){
    // sort edges by lexicografic non-descending endpoints
    // each edge has been created such that ep[0] < ep[1]
    std::sort(edges.begin(), edges.end(), CHAMedge::vertexSortFunc);

    for(uint32_t ei=0; ei < edges.size()-1; ){
        CHAMedge& e = edges[ei];

        while( ++ei < edges.size() && e.ep[0]==edges[ei].ep[0] && e.ep[1]==edges[ei].ep[1]){
            CHAMedge& ne = edges[ei];
            // each edge initially has only one incident triangle, ne never been merged
            e.inc_face.push_back( ne.inc_face[0] );
            ne.inc_face.clear(); // ne is now isolated
        }
    }

    // remove duplicated edges (no inc_faces) from edges vectcor
    edges.erase(std::remove_if(edges.begin(), edges.end(), CHAMedge::isIsolatedPtr), edges.end());
}

// Fill the PLCc data structure by using the input triangulation information
void PLCc::initialize(){
    
    // -- Verices -- 
    vertices.reserve(n_in_vrts);
    for(uint32_t i=0; i<n_in_vrts; i++) {
        const double* x = plc.coordinates.data() + i*3;
        vertices.push_back( new explicitPoint(*x, *(x +1), *(x +2)) ); 
    }
    mark_vrts.resize(vertices.size(),0);
    ref_exp3D_vrt.resize(vertices.size(), UINT32_MAX);

    // -- Edges --
    edges.resize(plc.numTriangles() * 3); // More than necessary, will be merged later

    // Fill only EV and partial ET relations
    // also fill vt (input_vertex - incident_input_triangle relation)
    input_vt.resize(vertices.size());
    uint32_t pei, v0,v1,v2;
    const uint32_t* t_ref = plc.triangle_vertices.data();
    for(uint32_t ti=0; ti<plc.numTriangles(); ti++){
        pei = ti*3;
        v0 = *(t_ref + pei   );  
        v1 = *(t_ref + pei +1);  
        v2 = *(t_ref + pei +2);
    
        // pre-edges
        edges[pei  ].fill_preEdge(v0, v1, ti);
        edges[pei+1].fill_preEdge(v1, v2, ti);
        edges[pei+2].fill_preEdge(v2, v0, ti);

        // vt
        input_vt[v0].push_back(ti);
        input_vt[v1].push_back(ti);
        input_vt[v2].push_back(ti);
    }

    // merge duplicated pre-edges and fill edges
    mergePreEdges(); // Now edges contains proper edges.
    mark_edges.resize(edges.size(), 0);

    // -- Faces --
    // Faces are input triangles

    // Find flat edges
    for(CHAMedge& e : edges){
        if( findIF_flat_edge(e) ) e.type = CHAMedge_t::flat; 
        else                      e.type = CHAMedge_t::undet;

        if (e.inc_face.size() & 1){ def_interior = false;}
    }

    #ifdef PLCC_VERBOSE_DEBUG
    disp_howManyFlatEdges();
    #endif

    // Create a face for each input triangle
    faces.resize( plc.numTriangles() );
    mark_faces.resize( faces.size(), 0);
    for(uint32_t i=0; i<plc.numTriangles(); i++) faces[i].triangle = i;

    // associate to each face (input triangle) its edges
    for(uint32_t ei=0; ei<edges.size(); ei++)
        for(uint32_t& t : edges[ei].inc_face) 
            faces[t].bounding_edges.push_back(ei);

    for(CHAMface& f : faces) orient_initial_triface_bnd(f);
}

// -------------------------- //
// SEARCH ACUTE PLCC ELEMENTS //
// -------------------------- //

// Assumes that the PLCedge e is one of the sides of the input-triangle ti.
// Returns the index of the vertex of ti different from e endpoints.
uint32_t PLCc::inTri_opp_vrt(const CHAMedge& e, const uint32_t ti) const {
  const uint32_t* tri_vrt = plc.triangle_vertices.data();
  uint32_t v = *( tri_vrt + (3*ti) ); if( !e.hasVertex(v) ) return v;
  v = *( tri_vrt + (3*ti +1) );       if( !e.hasVertex(v) ) return v;
  v = *( tri_vrt + (3*ti +2) );       if( !e.hasVertex(v) ) return v;
  return UINT32_MAX;
}

void PLCc::inTri_opp_edge(const uint32_t v, const uint32_t ti, uint32_t& u1, uint32_t& u2) const{
    if(v == plc.triangle_vertices[3*ti]){ 
        u1 = plc.triangle_vertices[3*ti +1];
        u2 = plc.triangle_vertices[3*ti +2];
        return;
    }

    if(v == plc.triangle_vertices[3*ti+1]){ 
        u1 = plc.triangle_vertices[3*ti +2];
        u2 = plc.triangle_vertices[3*ti ];
        return;
    }

    if(v == plc.triangle_vertices[3*ti+2]){ 
        u1 = plc.triangle_vertices[3*ti ];
        u2 = plc.triangle_vertices[3*ti +1];
        return;
    }

   assert(false); // SHOULD NEVER BEEN REACHED
}

// Returns true if the two incident faces at ei forms an acute dihedral angle.
bool PLCc::findIF_acute_edge(const CHAMedge& e) const {

    if(e.inc_face.size() == 1) return false;
    //if(e.inc_face.size() > 4) return true; // 369/5 < 360/4 = 90 
    if(e.inc_face.size() > 2) return true; // This line replaces the previous as the Delaunay Refinement 
                                           // algorithm (tetmesh.h) probably does not support non-manyfold configurations.

    uint32_t u, v;
    const pointType* e0_pt = vertices[ e.ep[0] ];
    const pointType* e1_pt = vertices[ e.ep[1] ];
    for(size_t i=1; i<e.inc_face.size(); i++){
        for(size_t j=0; j<i; j++){
            u = inTri_opp_vrt(e, faces[ e.inc_face[i] ].triangle);
            v = inTri_opp_vrt(e, faces[ e.inc_face[j] ].triangle);
            if( isAcuteDihedral_exact(e0_pt, e1_pt, vertices[u], vertices[v]) ) return true;
        }
    }

    return false;
}

// Returns true if at least:
// two non-flat edges incident at vi forms an acute angle, or
// an edge and an input triangle incident at vi forms an acute angle. 
bool PLCc::findIF_acute_vrt(const uint32_t vi, 
                            const std::vector<uint32_t>& vv_i, 
                            const std::vector<uint32_t>& vt_i    ) const {
    
    const pointType* vip = vertices[vi];

    // Check each edge-edge angles
    for(uint32_t j = 1; j < vv_i.size(); j++){
        for(uint32_t k = 0; k < j; k++){
            // isAcuteAngle is exact (see PLC.hpp), returns TRUE if acute at second argument
            if( isAcuteAngle( vertices[vv_i[j]], vip, vertices[vv_i[k]] ) ) return true;
        }
    }

    // Check angles formed by each edge and triangle sharing exactly a vertex.
    uint32_t t0, t1, t2, ep;
    for(uint32_t j = 0; j < vv_i.size(); j++){
        ep = vv_i[j];
        for(uint32_t k = 0; k < vt_i.size(); k++){
            uint32_t tid = vt_i[k] * 3;
            t0 = plc.triangle_vertices[tid   ];
            t1 = plc.triangle_vertices[tid +1];
            t2 = plc.triangle_vertices[tid +2];

            // vi is a vertex of vt_i[k], 
            // if ep is a vertex of vt_i[k] too
            // then <vi,ep> is an edge of vt_i[k]
            if( t0==ep || t1==ep || t2==ep ) continue;
            
            if( t1 == vi ) std::swap(t0,t1); 
            else if( t2 == vi ) std::swap(t0,t2); 
            // now t0 is vi

            if( isAcuteAngle( vertices[ep], vip, vertices[t1], vertices[t2] ) ) return true;
        }
    }    

  return false;
}

// returns the length of the shortest segment <vi, vv_i[j]>
double PLCc::closest_vv_dist(const uint32_t vi, const std::vector<uint32_t>& vv_i ) const {
    double sq_dist = vPtsSqDist(vi, vv_i[0]);
    double tmp_sq_dist;
    for(size_t j=1; j<vv_i.size(); j++){
        tmp_sq_dist = vPtsSqDist(vi, vv_i[j]); 
        if(tmp_sq_dist < sq_dist) sq_dist = tmp_sq_dist;
    }
    return sqrt( sq_dist );
}

// returns the shortest between the length of the segments that connect
// a vertex vi with the opposite sides of its incident triangles
double PLCc::closest_vt_dist(const uint32_t vi, const std::vector<uint32_t>& vt_i ) const {
    uint32_t u1, u2;   inTri_opp_edge(vi, vt_i[0], u1, u2);
    double sq_dist = vSqVrtDistSeg(vi, u1, u2);
    double tmp_sq_dist;
    for(size_t j=0; j<vt_i.size(); j++){
        inTri_opp_edge(vi, vt_i[j], u1, u2);
        tmp_sq_dist = vSqVrtDistSeg(vi,u1,u2);
        if( tmp_sq_dist < sq_dist ) sq_dist = tmp_sq_dist;
    }
    return sqrt( sq_dist );
}

// Search for acute vertices and edges and computes ch_dist.
// For acute vertices vrt_ch_dist is 1/3 of the shortest incident edge length.
void PLCc::search_acute_angles(){

    if( vrt_ch_dist.size() != vertices.size() ) vrt_ch_dist.resize(vertices.size(), -1.0); // init

    // Classify not-flat edges between acute and not-acute.
    std::vector< std::vector<uint32_t> > vv(vertices.size()); // v -> all v' connected with v by a not-flat edge
    for(CHAMedge& e : edges) if( !e.isFlat() ){
        vv[e.ep[0]].push_back(e.ep[1]);
        vv[e.ep[1]].push_back(e.ep[0]);

        if( findIF_acute_edge(e) ) e.type = CHAMedge_t::acute;
    }

    #ifdef PLCC_VERBOSE_DEBUG
    disp_howManyAcuteEdges();
    #endif
    #ifdef PLCC_VERBOSE_DEBUG_LEV1
    disp_allAcuteEdges();
    #endif

    // Force acutness of acute edge endpoints (necessary for edge chamfering)
    for(const CHAMedge& e : edges) if(e.isAcute()){ 
        uint32_t e0 = e.ep[0], e1 = e.ep[1];
        if(!is_acute_vrt(e0)) vrt_ch_dist[e0] = get_vrt_ch_dist( e0, vv[e0], input_vt[e0] );
        if(!is_acute_vrt(e1)) vrt_ch_dist[e1] = get_vrt_ch_dist( e1, vv[e1], input_vt[e1] );
    }

    // Classify acute vertices
    for(uint32_t vi=0; vi<vertices.size(); vi++) if( !is_acute_vrt(vi) ) {
        if( findIF_acute_vrt(vi, vv[vi], input_vt[vi]) ) 
            vrt_ch_dist[vi] = get_vrt_ch_dist( vi, vv[vi], input_vt[vi] );
    }

    // Chamfering distance must take in to account flat edge lenght,
    // otherwise flat edges may be "cutted out". 
    for(const CHAMedge& e : edges ) if( e.isFlat() ){
        uint32_t e0 = e.ep[0], e1 = e.ep[1];
        double e_len = sqrt( vPtsSqDist(e0, e1) );
        if( vrt_ch_dist[e0] * 3 > e_len ) vrt_ch_dist[e0] = e_len / 3;
        if( vrt_ch_dist[e1] * 3 > e_len ) vrt_ch_dist[e1] = e_len / 3;
    }

    #ifdef PLCC_VERBOSE_DEBUG
    disp_howManyAcuteVertices();
    #endif
    #ifdef PLCC_VERBOSE_DEBUG_LEV1
    disp_allAcuteVertices();
    #endif
}

// ----------------------------------------------------- //
// CHAMFERING (additional) POINTS CREATION AND PLACEMENT //
// ----------------------------------------------------- //

// New BPT has as reference-triangle <p0, s0, s1> where s0 and s1 are the original vertices
// associated to the LNC p1 (i.e. p1 = s0 * (1-tau) + si * tau)
// implicitPoint3D_BPT* get_BPTonSegm_with_ExpLNC_ep(const explicitPoint3D& p0, const implicitPoint3D_LNC& p1, const double t){
//     const explicitPoint3D& s0 = p1.P();
//     const explicitPoint3D& s1 = p1.Q();
//     const double tau = p1.T();
//     return new implicitPoint3D_BPT(s0, s1, p0, t*(1-tau), t*tau);
// }

// ASSUMPTIONS:
// p0 = s0 * (1-teta) + s2 * teta
// p1 = s0 * (1-tau ) + s1 * tau
// i.e. LNC points share a common parent vertex s0 which is stored in the field P() of both of them.
// New BPT has as reference-triangle <s0, s1, s2>
implicitPoint3D_BPT* get_BPTonSegm_with_LNCLNC_ep_baseCase(const implicitPoint3D_LNC& p0, const implicitPoint3D_LNC& p1, const double t){
    const explicitPoint3D& s0 = p0.P();
    const explicitPoint3D& s1 = p1.Q();
    const explicitPoint3D& s2 = p0.Q();
    const double teta = p0.T();
    const double tau = p1.T();
    return new implicitPoint3D_BPT(s1, s2, s0, t*tau, (1-t)*teta);
}

// ASSUMPTION: one of the reference-points of p0 and p1 must be shared between them.
// New BPT has as reference-triangle <s0, s1, s2>
implicitPoint3D_BPT* get_BPTonSegm_with_LNCLNC_ep(const implicitPoint3D_LNC& p0, const implicitPoint3D_LNC& p1, const double t){

    if( p0.P() == p1.P() ) return get_BPTonSegm_with_LNCLNC_ep_baseCase(p0, p1, t);
    else if( p0.P() == p1.Q() ){
        const implicitPoint3D_LNC swap_p1(p1.Q(), p1.P(), 1-(p1.T()) );
        return get_BPTonSegm_with_LNCLNC_ep_baseCase(p0, swap_p1, t);
    }
    else if( p0.Q() == p1.P() ){
        const implicitPoint3D_LNC swap_p0(p0.Q(), p0.P(), 1-(p0.T()));
        return get_BPTonSegm_with_LNCLNC_ep_baseCase(swap_p0, p1, t);
    }
        
    assert( p0.Q() == p1.Q() );
    const implicitPoint3D_LNC swap_p0(p0.Q(), p0.P(), 1-(p0.T()));
    const implicitPoint3D_LNC swap_p1(p1.Q(), p1.P(), 1-(p1.T()));
    return get_BPTonSegm_with_LNCLNC_ep_baseCase(swap_p0, swap_p1, t);
}

// Creates a new point p on the segment <v0,v1> at distance d from:
// - v0 if d0 = 0
// - v1 if d0 = 1
// The point is EXACTLY on the segment (implicit point) when implicit point are used (no macro FP_CHAMFERING)
uint32_t PLCc::new_vrt_on_segment(uint32_t v0, uint32_t v1, const double d, const uint32_t d0){
    
    if(d0==1) std::swap(v0, v1);
    double l = sqrt( vPtsSqDist(v0, v1) );
    double t = d / l;

    #ifdef PLCC_DEBUG
    if( t < 0 || t > 1 ) std::cout<<"[cham.cpp new_vrt_on_segment()] ERROR: d = "<<d<<" - l = "<<l<<"\n";
    assert( t > 0 && t < 1 );
    #else
    if( t <= 0 || t >= 1 ) std::cout<<"[cham.cpp new_vrt_on_segment()] WARNING: invalid t = "<<t<<" for LNC; d = "<<d<<" - l = "<<l<<"\n";
    #endif
  
    // new point p is such that p = (1-t) v0 + t v1
    pointType* p;
    uint32_t exp3D_i = (t < 0.5) ? v0 : v1;

    #ifndef FP_CHAMFERING

    const pointType* p0 = vertices[v0];
    const pointType* p1 = vertices[v1];

    if(p0->isExplicit3D() && p1->isExplicit3D()){
        p = new implicitPoint_LNC( p0->toExplicit3D(), p1->toExplicit3D(), t);
    }
    else if(p0->isLNC() && p1->isLNC()){
        std::cout<<"[cham.cpp - new_vrt_on_segment()] WARNING should pass here only using fallbaxk solution\n";
        p = get_BPTonSegm_with_LNCLNC_ep( p0->toLNC(), p1->toLNC(), t);
    }
    else{
        // THIS SHOLD NEVER BE REACHED
        #ifdef TEST_CHAMFERING
        exit( (int) EXIT_t::vons_fail );
        #endif
        #ifdef PLCC_DEBUG
        if(d0==1) std::swap(v0, v1);
        report_vons_fail(v0, v1);
        assert( false );
        #endif
        // ip_error("[cham.cpp - new_vrt_on_segment()] ERROR invalid edge endpoints\n");
        std::cout<<"[cham.cpp - new_vrt_on_segment()] ERROR Invalid edge endpoints\n"; exit(1);
    }

    #else

    // NOT implictit -> use FP approximation
    // const vector3d Ov0( vertices[v0] );
    // const vector3d Ov1( vertices[v1] );
    // const vector3d Op( Ov0 * (1-t) + Ov1 * t );
    const vector3d Op = vector3d( vertices[v0] ).leftLinComb( vector3d( vertices[v1] ), t );
    p = new explicitPoint3D(Op.c[0], Op.c[1], Op.c[2]);

    #endif

    add_vertex( p, exp3D_i);
    return ((uint32_t) vertices.size()-1);
}

// DEBUG
void display_anomalous_bpt_coeff_data(
                double& xi0, double& eta0, double& xi1, double& eta1, 
                const vector3d& Ov, const vector3d& Ou0, const vector3d& Ou1, 
                double t0, double t1,
                double d, const double zero_toll,
                double cos_alpha, 
                double teta, double psi, double psi_lim){

    vector3d Op0 = Ov.leftLinComb(Ou0,t0);
    vector3d Op1 = Ov.leftLinComb(Ou1,t1);
    vector3d Oq0 = (Ou0 - Ov) * xi0 + (Ou1 - Ov) * eta0 + Ov;
    vector3d Oq1 = (Ou0 - Ov) * xi1 + (Ou1 - Ov) * eta1 + Ov;
    double dist = (Oq0-Oq1).sq_length();
    if(dist < zero_toll*zero_toll){
        std::cout<<"d="<<d<<", cos_alpha="<<cos_alpha<<"\n"; 
        std::cout<<"||q1-q0|2 - |p0-q0|2|="<<abs(dist - (Op0 - Oq0).sq_length())<<"\n";
        std::cout<<"||q1-q0|2 - |p1-q1|2|="<<abs(dist - (Op1 - Oq1).sq_length())<<"\n";
        std::cout<<"||p0-q0|2 - |p1-q1|2|="<<abs((Op0 - Oq0).sq_length() - (Op1 - Oq1).sq_length())<<"\n";
        std::cout<<"|p0-v|2="<<(Op0 - Ov).sq_length()<<"\n";
        std::cout<<"|p1-v|2="<<(Op1 - Ov).sq_length()<<"\n";
        std::cout<<"v = "<<Ov<<"\n";
        std::cout<<"t0 = "<<t0<<", u0 = "<<Ou0<<"\n";
        std::cout<<"t1 = "<<t1<<", u1 = "<<Ou1<<"\n";
        std::cout<<"q0 = "<<Oq0<<"\n";
        std::cout<<"q1 = "<<Oq1<<"\n";
        std::cout<<"xi0 = "<<xi0<<", eta0 = "<<eta0<<"\n";
        std::cout<<"xi1 = "<<xi1<<", eta1 = "<<eta1<<"\n";
        std::cout<<"psi = "<<psi<<"\n";
        std::cout<<"psi_lim = "<<psi_lim<<"\n";
    }
}

// Computes the cosine of the angle at A of the triangle <A,B,C>
// using the cosine law.
double tri_oppA_cos(const vector3d& A, const vector3d& B, const vector3d& C){
    double a_sq = B.dist_sq( C );
    double b_sq = C.dist_sq( A );
    double c_sq = A.dist_sq( B );
    return (b_sq + c_sq - a_sq) / (2 * sqrt(b_sq * c_sq) );
}

// Computes the baricentric coordinates (xi0,eta0), (xi1,eta1) of
// two BPT to replace segments <LNC0,V> <V,LNC1> with a "bridge"
// <LNC0,BPT0>, <BPT0,BPT1>, <BPT1,LNC2> in orther that
// <LNCi,BPTi> is orthognal to <LNCi,V> for i=0,1.
// The function try to equalize new segments, i.e.
// dist(LNC1,BPT1) = dist(BPT1,BPT2) = dist(BPT2,LNC2)
bool get_bpt_coeff(double& xi0, double& eta0, double& xi1, double& eta1, 
                   const vector3d& Ov, 
                   const vector3d& Ou0, const vector3d& Ou1, 
                   double t0, double t1,
                   double d, const double zero_toll = 0.000000000000001){

    double cos_alpha = tri_oppA_cos(Ov, Ou0, Ou1); // cos at v

    double one_plus_cos_alpha = 1+cos_alpha;
    double one_min_cos_alpha = 1-cos_alpha;

    if(abs( one_plus_cos_alpha ) < zero_toll ||
       abs( one_min_cos_alpha ) < zero_toll    ){ 
        std::cout<<"WARNING: [get_bpt_coeff] parameters below zero_toll\n";
        //return false;
    }

    double teta = 1.0 / one_plus_cos_alpha;
    double psi, psi_epsilon;

    // optimal psi: generated BPT produces 3 equaly long segments with
    // k = ( d * sqrt(2) * sqrt(1 - cos_alpha) )  /  ( 1 + sqrt(2) * sqrt(1 + cos_alpha) ).
    // psi = (k / d) * sqrt(  (1 + cos_alpha) / (1 - cos_alpha ) ), thus...
    psi = 1.0  /  ( 1.0 + 1.0 / sqrt(one_plus_cos_alpha * 2.0) );
    psi_epsilon = sqrt( one_plus_cos_alpha / one_min_cos_alpha ); // this is psi when k = d
    if(psi > psi_epsilon) psi = psi_epsilon; // BTPs have to stay at most at distance d from triangle side.

    // if teta * (t0+t1) < 1.0 THEN BTPs SHURELY stay inside the triangle.
    double psi_lim = -1.0;
    if( teta * (t0+t1) > 1.0 ){ 
        // BPTs may stay otside triangles
        double psi_lim_0 = (1.0-t0) / ( (teta-1.0)*t0 + teta*t1 );
        double psi_lim_1 = (1.0-t1) / ( (teta-1.0)*t1 + teta*t0 );
        psi_lim = 0.5 * min(psi_lim_0, psi_lim_1);
        if( psi > psi_lim ) psi = psi_lim; // at least one BPT is outside the triangle
    }
    
    double psi_teta = psi * teta;
    double term = 1.0 + psi * (teta - 1.0);

    xi0 = t0 * term;
    eta0 = t1 * psi_teta;
    xi1 = t0 * psi_teta;
    eta1 = t1 * term;

    // DEBUG
    display_anomalous_bpt_coeff_data(xi0,eta0,xi1,eta1, Ov,Ou0,Ou1, t0,t1, d,zero_toll,cos_alpha, teta, psi, psi_lim);

    bool not_representable = (
        xi0  <= 0.0 || xi0  >= 1.0 ||
        eta0 <= 0.0 || eta0 >= 1.0 ||
        xi1  <= 0.0 || xi1  >= 1.0 ||
        eta1 <= 0.0 || eta1 >= 1.0 ||
        xi0 + eta0 >=1.0 || xi1 + eta1 >=1.0 );

    if( not_representable ){ 

        // DEBUG
        std::cout<<"WARNING: [get_bpt_coeff] Steiner point not representable as BPT\n";
        std::cout<<"xi0 = "<<xi0<<", eta0 = "<<eta0<<"\n";
        std::cout<<"xi1 = "<<xi1<<", eta1 = "<<eta1<<"\n";
        std::cout<<"d = "<<d<<", k = "<<psi*d*sqrt(one_min_cos_alpha/one_plus_cos_alpha)<<", t0 = "<<t0<<", t1 = "<<t1<<"\n"; 
        std::cout<<"cos_alpha = "<<cos_alpha<<"\n";
        std::cout<<"teta = "<<teta<<", psi = "<<psi<<", psi_lim = "<<psi_lim<<"\n";

        return false; 
    }

    return true;
}

// Creates a 2 new points p and p' inside the input triangle indexed as itri_ind.
// The point is EXACTLY inside the triangle (implicit point) when macro FP_CHAMFERING is commented,
// otherwise floating point coordinates are used to represent p.
uint32_t PLCc::new_vrts_in_inputTri(const uint32_t fi, 
                                    const uint32_t vi, 
                                    const uint32_t u0, const uint32_t ou0, 
                                    const uint32_t u1, const uint32_t ou1){

    #ifdef PLCC_DEBUG
    assert( vertices[ou0]->isExplicit3D() );
    assert( vertices[ou1]->isExplicit3D() );
    #endif

    const vector3d Ov( vertices[vi] );
    const vector3d Ou0( vertices[ou0] );
    const vector3d Ou1( vertices[ou1] );

    double xi0, eta0, xi1, eta1, t0, t1;

    #ifndef FP_CHAMFERING

    assert(vertices[u0]->isLNC() && vertices[u1]->isLNC());

    t0 = vertices[u0]->toLNC().T();
    t1 = vertices[u1]->toLNC().T();

    #else
    
    t0 = sqrt( (Ov - vector3d(vertices[u0])).sq_length() / (Ov - vector3d(vertices[ou0])).sq_length() );
    t1 = sqrt( (Ov - vector3d(vertices[u1])).sq_length() / (Ov - vector3d(vertices[ou1])).sq_length() );
    
    #endif

    double r = vrt_ch_dist[vi]; if(epsilon < r) r = epsilon;
    if(r<0) r = epsilon; // ch_dist may be -1 if not initialized
    if( !get_bpt_coeff(xi0,eta0, xi1,eta1, Ov,Ou0,Ou1, t0,t1, r) ) return INVALID_BPT;
  
    pointType *p0, *p1;
    
    #ifndef FP_CHAMFERING

    const explicitPoint3D& pv = vertices[vi]->toExplicit3D();
    const explicitPoint3D& pu0 = vertices[ou0]->toExplicit3D();
    const explicitPoint3D& pu1 = vertices[ou1]->toExplicit3D();

    p0 = new implicitPoint3D_BPT(pu0, pu1, pv, xi0, eta0);
    p1 = new implicitPoint3D_BPT(pu0, pu1, pv, xi1, eta1);

    #else

    // FP approximation
    const vector3d Op0( Ou0 * xi0 + Ou1 * eta0 + Ov * ( 1 - (xi0+eta0) ) );
    const vector3d Op1( Ou0 * xi1 + Ou1 * eta1 + Ov * ( 1 - (xi1+eta1) ) );
    p0 = new explicitPoint3D(Op0.c[0], Op0.c[1], Op0.c[2]);
    p1 = new explicitPoint3D(Op1.c[0], Op1.c[1], Op1.c[2]);

    #endif
    
    add_vertex( p0, vi );
    add_vertex( p1, vi );
    return ((uint32_t) vertices.size()-1);
}

// splt edge ei at distance d from endpoint edeges[ei].ep[ep_i]
void PLCc::chamfer_edge_ep(const size_t ei, double d, const uint32_t ep_i){

    CHAMedge& e = edges[ei]; // e = <e0,e1>
    const uint32_t e0 = e.ep[0];

    #ifdef PLCC_DEBUG
    std::cout<<"\nchamfering (d="<<d<<") ep "<<ep_i<<" of "; print_edge(ei);
    #endif

    const uint32_t Pt_i = new_vrt_on_segment(e.oep[0], e.oep[1], d, ep_i);
    
    // Update PLCedges
    CHAMedge new_e(e); // copy e
    e.ep[1] = Pt_i; // Update PLCedge e endpoints: <e0,e1> becomes <e0,Pt>
    new_e.ep[0] = Pt_i; // New edge is <Pt_i,e1>
    edges.push_back( new_e ); // <Pt,e1> 
    mark_edges.push_back(0);

    // Update type
    if(ep_i==1) edges.back().type = CHAMedge_t::junk;
    else        edges[ei].type = CHAMedge_t::junk;

    #ifdef PLCC_DEBUG
    std::cout<<"after chamfering:\n"; print_edge(ei); print_edge(edges.size()-1);
    if(ep_i==1){ std::cout<<"edge["<< (uint32_t) (edges.size()-1) <<"] is junk\n"; }
    else{ std::cout<<"edge["<< (uint32_t)ei <<"] is junk\n"; }
    #endif

    // update connectivity (boundaries of incident faces)
    uint32_t new_ei = (uint32_t) edges.size()-1;
    for(uint32_t fi : edges[ei].inc_face){

        faces[fi].make_last((uint32_t)ei); // <e0,Pt> is the last on face boundary.
        if(edges[ faces[fi].bounding_edges.front() ].hasVertex(e0)){
            faces[fi].make_first((uint32_t)ei);
        }
        faces[fi].bounding_edges.push_back( new_ei );

        #ifdef PLCC_DEBUG
        std::cout<<"\n"; print_face_edges(fi);
        #endif
    }
}

// Splits each edge that has at least one acute vertex v as endpoint,
// Steiner point is inserted at distance r = min(epsilon, lfs(v)/3) from v.
// NOTE. Steiner points have only the two sub-edges as incident edges
//       thus they are not acute.
void PLCc::chamfering_vrts(){
    size_t n_orig_edges = edges.size();
    uint32_t v;  double r;
    for(size_t ei=0; ei<n_orig_edges; ei++) if(!edges[ei].isIsolated()){

        v = edges[ ei ].ep[ 1 ];
        if( is_acute_vrt(v) ){ 
            r = vrt_ch_dist[v]; if(epsilon < r) r = epsilon;
            chamfer_edge_ep(ei, r, 1); // <- diveides edges[ei] at distance r from endpoint 1
        }
        // now edges[ei] is <e0,e1> or <e0,new_point>

        v = edges[ ei ].ep[ 0 ];
        if( is_acute_vrt(v) ){ 
            r = vrt_ch_dist[v]; if(epsilon < r) r = epsilon;
            chamfer_edge_ep(ei, r, 0);
        }
    }
}

// ---------------- //
// BASIC CHAMFERING //
// ---------------- //

void PLCc::remove_junk_edges_from_face(uint32_t fi){

    #ifdef PLCC_DEBUG
    std::cout<<"before removing..\n"; print_face_edges(fi);
    #endif

    CHAMface& f = faces[fi];
    const std::vector<uint32_t>& fbnd = f.bounding_edges;
    // if a face is bounded by a junk edge,
    // there must be another consecutive junk edge
    // that share with the first an acute vertex 

    for(size_t i=0; i<fbnd.size(); i++){
        uint32_t ei = fbnd[i];
        if( ei == UINT32_MAX || !edges[ ei ].isJunk() ) continue;
        
        const CHAMedge& e1 = edges[ ei ]; // is a junk edge
        
        // 1) find its consecutive edge sharing an acute vertex
        size_t j, prec=i, next=i; 
        f.advance_on_bnd(next);
        if( fbnd[next] != UINT32_MAX || 
            !edges[ fbnd[next] ].isJunk() || 
            !e1.has_commonVertex( edges[ fbnd[next] ] ) ){
            j = next;
        }
        else if( fbnd[prec] != UINT32_MAX || 
                 !edges[ fbnd[prec] ].isJunk() || 
                 !e1.has_commonVertex( edges[ fbnd[prec] ] ) ){ 
            j = prec;
        }

        #ifdef PLCC_DEBUG
        assert(j != i);   
        #endif

        // 2) create a new edge connecting non-acute vertices of e1 and e2
        const CHAMedge& e2 = edges[ fbnd[j] ];
        uint32_t v1 = e1.notCommonVertex(e2);
        uint32_t v2 = e2.notCommonVertex(e1);
        edges.push_back( CHAMedge(v1, v2, fi) );
        mark_edges.push_back(0);

        // 3) Replace e1 with e_new and remove e2 from f
        f.bounding_edges[i] = (uint32_t) edges.size() -1;
        f.bounding_edges[j] = UINT32_MAX;
    }

    f.bounding_edges.erase(
        std::remove_if( f.bounding_edges.begin(), 
                        f.bounding_edges.end(), 
                        [](uint32_t x){ return x==UINT32_MAX;}), 
        f.bounding_edges.end());

    #ifdef PLCC_DEBUG
    std::cout<<"after removing..\n"; print_face_edges(fi);
    #endif
}

//
void get_backProjection_dist(const pointType* el, const pointType* er, 
                             const pointType* ul, const pointType* ur, 
                             double d, double& dl, double& dr){
    
    // Given:
    // - the segments <ul,el>, <el,er> and <er,ur> such that:
    //   * ul, el, er, ur are coplanar;
    //   * ul, ur belong to the same half-plane defined by the straight line for el and er;
    //   * the angles ul - el - er and el - er - ur are greater than pi/2;
    // - the distance d.
    // Consider the stright line s parallel to <el,er> and distant d from <el,er>.
    // intersecting <ul,el> and <er,ur> in pl and pr respectively.
    // We want to detrmine:
    // - dl = dist( el, pl )
    // - dr = dist( er, pr ) 

    // Let be:
    //    alpha = angle(ul - el - er) - pi/2 [alpha is acute],
    //    beta  = angle(el - er - ur) - pi/2 [beta  is acute],
    //    l_vett = ul - el
    //    c_vett = er - el
    //    c_ovett = el - er
    //    r_vett = ur - er

    // Following relations holds (for <ul,el> but the same hold when <er,ur> is considered) :
    //
    // dl = d / cos( alpha )
    //
    // cos( alpha + pi/2 ) = < l_vett * c_vett > / ( norm(l_vett) * norm(c_vett) )
    //
    // sin^2( alpha ) = < l_vett * c_vett >^2 / ( norm^2(l_vett) * norm^2(c_vett) )
    // cos^2( alpha ) = 1 - < l_vett * c_vett >^2 / ( norm^2(l_vett) * norm^2(c_vett) )
    // norm^2(l_vett) * norm^2(c_vett) * cos^2( alpha ) = norm^2(l_vett) * norm^2(c_vett) - < l_vett * c_vett >^2
    //
    // define: den_l = norm^2(l_vett) * norm^2(c_vett) - < l_vett * c_vett >^2
    //
    // dl^2 = d^2 * norm^2(l_vett) * norm^2(c_vett) / den_l


    // NOTE. HYPOTESIS (for <ul,el> but the same hold when <er,ur> is considered) 
    // 0 < alpha < pi/2 => 0 < cos( alpha ) < 1 
    // thus den_l != 0 

    // NOTE. STABILITY (for <ul,el> but the same hold when <er,ur> is considered)
    // dr or dl computation is unstable when alpha is close to -pi/2 or pi/2 as cos(aplha) goes to 0.
    // 1) alpha > 0 thus alpha >> -pi/2; (NO PROBLEM) 
    // 2) UNFORTUNATELY alpha may be arbitraly close to pi/2; (PROBLEM)

    // SOLUTION:
    // In that cases (STABILITY - 2) dl is not well defined, but we can
    // take as dl 1/3 * length(<ul,el>) and try to cmpute a choerent dr,
    // unfortunately being new_d = 1/3 * length(<ul,el>) * sin(alpha) very small too
    // it will be probabily necessary to compute dr with the original d.
    // As consequence <pl,pr> will not be parallel to <el,er>:
    // the angles ul - pl - pr and pl - pr -ur will be different from alpha and beta
    // but they will be greater of alpha and beta respectivelly.

    dl = -1.0;
    dr = -1.0;

    const vector3d Oepr( er );
    const vector3d Oepl( el );
    const vector3d e_vect( Oepl - Oepr );
    const vector3d er_vect( vector3d( ur ) - Oepr );
    const vector3d el_vect( vector3d( ul ) - Oepl );
    double sq_len_e = e_vect.sq_length();
    double sq_len_er = er_vect.sq_length();
    double sq_len_el = el_vect.sq_length();
    double norm_prod_r = sq_len_e * sq_len_er;
    double norm_prod_l = sq_len_e * sq_len_el;
    double e_dot_er = e_vect.dot(er_vect);
    double e_dot_el = e_vect.dot(el_vect); // it would be -e_vect dot el_vect but then is squared so dont mind..

    double den_r = norm_prod_r - e_dot_er*e_dot_er;
    double den_l = norm_prod_l - e_dot_el*e_dot_el;
    double r_coeff_sq, l_coeff_sq;

    double toll = 0.001;
    if( den_r < toll ) dr = sqrt( sq_len_er )/3;
    else r_coeff_sq = norm_prod_r / den_r;

    if( den_l < toll ) dl = sqrt( sq_len_el )/3;
    else l_coeff_sq = norm_prod_l / den_l;

    if( dr != -1.0 && dl != -1.0 ) return;

    double d_sq = d*d;
    double dl_sq, dr_sq = d_sq * r_coeff_sq;

    if( dl == -1.0 ){ 
        dl_sq = d_sq * l_coeff_sq;
        if( dl_sq > sq_len_el/9 ) dl = sqrt( sq_len_el )/3;
    }

    if( dr == -1.0 ){ 
        dr_sq = d_sq * r_coeff_sq;
        if( dr_sq > sq_len_er/9 ) dr = sqrt( sq_len_er )/3;
    }

    if( dl == -1.0 ) dl = d * sqrt(l_coeff_sq);
    if( dr == -1.0 ) dr = d * sqrt(r_coeff_sq);

    // PARALLELISM MAY BE CONSERVED IN MORE SITUATIONS..
    // double d_sq_rlim = -1.0, d_sq_llim = -1.0;

    // if( dr == -1.0 && sq_len_er/9 < dr_sq ){
    //     d_sq_rlim = (norm_prod_r - e_dot_er*e_dot_er)/(9 * sq_len_e);
    // }

    // if( dr == -1.0 && sq_len_el/9 < d_sq*l_coeff_sq ){
    //     double d_sq_llim = (sq_len_el/9) * (1 - (e_dot_er*e_dot_er)/norm_prod_l );
    //     dl = sqrt( sq_len_el )/3;
    //     return;
    // }

    // if( sq_len_el/9 < d_sq*l_coeff_sq || sq_len_er/9 < d_sq*r_coeff_sq ){
    //     double d_sq_llim = (sq_len_er/9) * (1 - (e_dot_el*e_dot_el)/norm_prod_l );
    //     double d_sq_rlim = (sq_len_er/9) * (1 - (e_dot_er*e_dot_er)/norm_prod_r );
    //     if(d_sq_llim < d_sq_rlim){
    //         dr = sqrt( d_sq_llim * r_coeff_sq );
    //         dl = sqrt( sq_len_er )/3;
    //     }
    //     else{
    //         dr = sqrt( sq_len_er )/3;
    //         dl = sqrt( d_sq_rlim * l_coeff_sq );
    //     } 
    // }
    // else{
    //     dr = d * sqrt( r_coeff_sq );
    //     dl = d * sqrt( l_coeff_sq );
    // }
}

//
void PLCc::chamfer_acute_edge_from_inc_face(uint32_t ei, uint32_t fi, double d){

    CHAMedge& e = edges[ ei ];
    CHAMface& f = faces[ fi ];

    #ifdef PLCC_DEBUG
    std::cout<<"\nacute "; print_edge(ei);
    std::cout<<"with "<<e.inc_face.size()<<" incident faces\n";
    std::cout<<"inc_face #"<<fi<<" -> "; print_face_edges(fi); 
    #endif

    // 1) Find er and el edges adjacent to e on boundary of f,
    //    find epr and epl the common endpoint between e-er and e-el respectivelly.
    //    er = <epr,ur> and el = <epl,ul>

    uint32_t eri, eli; // indices wrt edges of er and el
    u32vect_iter er_it, el_it; // iterators to er and el on boundary of f
    uint32_t epl, epr; // endpoints of er and el shared with e (or e_chain)
    uint32_t ul, ur; // endpoints of er and el not-shared with e (or e_chain)

    epl = e.ep[0], epr = e.ep[1];
    el_it = std::find(f.bounding_edges.begin(), f.bounding_edges.end(), ei);
    er_it = el_it; 

    u32vect_iter it=er_it; f.advance_on_bnd(it);
    if(e.commonVertex(edges[*it])==epl) std::swap(epl, epr);

    f.reverse_on_bnd(el_it);
    eli = *el_it;
    ul = edges[eli].oppositeVertex(epl);

    f.advance_on_bnd(er_it);
    eri = *er_it;
    ur = edges[eri].oppositeVertex(epr);

    CHAMedge& er = edges[eri];
    CHAMedge& el = edges[eli];

    #ifdef PLCC_DEBUG
    if(er.inc_face.size()!=1){
        std::cout<<"[chamfer_acute_edge_from_inc_face] ERROR: "
                 <<"edge "<<eri<<" has "<<er.inc_face.size()<<" inc_faces.\n";
        print_edge_and_inc_face(eri);
    }
    if(el.inc_face.size()!=1){
        std::cout<<"[chamfer_acute_edge_from_inc_face] ERROR: "
                 <<"edge "<<eli<<" has "<<el.inc_face.size()<<" inc_faces.\n";
        print_edge_and_inc_face(eli);
    }
    #endif

    assert(er.inc_face.size()==1);
    assert(el.inc_face.size()==1);

    double dr, dl;
    get_backProjection_dist(vertices[epl], vertices[epr], vertices[ul], vertices[ur], d, dl, dr);
   
    #ifdef PLCC_DEBUG
    std::cout<<" d = "<<d<<"; dl = "<<dl<<"; dr = "<<dr<<"\n";
    #endif

    // 3) Create new vertices vr and vl on er at distance dr from epr
    //    and on el at distance dl from epl respectivelly
    uint32_t vr_i, vl_i;
    uint32_t our = ur, oul = ul;
    if( vertices[ur]->isBPT() ) our = er.relatedOriginalVertex(ur);
    if( vertices[ul]->isBPT() ) oul = el.relatedOriginalVertex(ul);

    vr_i = new_vrt_on_segment(epr, our, dr, 0);
    vl_i = new_vrt_on_segment(epl, oul, dl, 0);

    // 4) Change endpoint epr of er with vr, change enpoint epl with vl.

    if( er.ep[0] == epr ) er.ep[0] = vr_i;
    else                  er.ep[1] = vr_i;

    if( el.ep[0] == epl ) el.ep[0] = vl_i;
    else                  el.ep[1] = vl_i;

    // 5) create a new edge en of endpoints vr and vl,
    //    replace e with en on boundary of f

    if(e.inc_face.size() > 1){
        edges.push_back( CHAMedge(vr_i, vl_i, fi) );
        mark_edges.push_back(0);
        f.replaceEdge_11(ei, (uint32_t) edges.size() -1);
        e.removeIncidentFace(fi);

        #ifdef PLCC_DEBUG
        std::cout<<"replaced by new "; print_edge(edges.size()-1);
        std::cout<<endl;
        #endif
    }
    else{
        e.ep[0] = vl_i;
        e.ep[1] = vr_i;
        e.oep[0] = e.ep[0];
        e.oep[1] = e.ep[1];
        e.type = CHAMedge_t::undet;
        assert( e.inc_face.size()==1 && e.inc_face[0]==fi );

        #ifdef PLCC_DEBUG
        std::cout<<"updated as "; print_edge(ei);
        std::cout<<endl;
        #endif
    }
    
}

//
// void PLCc::get_edge_ch_dist(std::vector<double>& edge_ch_dist) {
//     double sq_epsilon = epsilon * epsilon;
//     for(size_t i=0; i<edges.size(); i++) if(edges[i].isAcute()){
//         const CHAMedge& e = edges[i];
//         double edge_sqlen = vPtsSqDist(e.ep[0], e.ep[1]);
//         edge_ch_dist[i] = min( sq_epsilon, edge_sqlen );
//         for(uint32_t fi : e.inc_face){
//             std::vector<uint32_t> fv; get_face_vertices( faces[fi], fv );
//             for(uint32_t vi : fv) if( !e.hasVertex(vi) )
//                 edge_ch_dist[i] = min( edge_ch_dist[i], vSqVrtDistEdge(vi, e) );
//         }
//         edge_ch_dist[i] = sqrt(edge_ch_dist[i])/3;
//         assert(edge_ch_dist[i] > 0);
//     } 
// }

// Replace face fi which is assumed to be contained in or coincident to
// its related input triangle ti (faces[fi].triangle) with
// a "basic" chamfered version of ti.
void PLCc::inputTriangleChamfering(uint32_t fi){
    const CHAMface& f = faces[fi];

    // Isolate f from incident faces:
    // 1- remove fi from all edges bounding f
    for(const uint32_t& ei : f.bounding_edges) if(ei!=UINT32_MAX){ 
        edges[ ei ].removeIncidentFace(fi);
    }

    // 2- substituite f with its related input triangle
    uint32_t ti = f.triangle;
    const uint32_t* in_tv = plc.triangle_vertices.data() + ti*3;
    uint32_t v[] = { *in_tv, *(in_tv +1), *(in_tv +2) };
    edges.push_back( CHAMedge(v[0], v[1], fi ) );
    edges.push_back( CHAMedge(v[1], v[2], fi ) );
    edges.push_back( CHAMedge(v[2], v[0], fi ) );
    mark_edges.push_back(0);
    mark_edges.push_back(0);
    mark_edges.push_back(0);
    
    faces[fi].bounding_edges.clear();
    uint32_t l = (uint32_t) edges.size();
    uint32_t be[] = {--l, --l, --l};
    faces[fi].bounding_edges.assign( be, be+3 );
    
    orient_initial_triface_bnd(faces[fi]);

    // chamfer new face
    edges[be[0]].type = CHAMedge_t::acute;
    edges[be[1]].type = CHAMedge_t::acute;
    edges[be[2]].type = CHAMedge_t::acute;

    double lv0v1, lv1v2, lv2v0, dv0, dv1, dv2, d;
    vector3d Ov0( vertices[v[0]] );
    vector3d Ov1( vertices[v[1]] );
    vector3d Ov2( vertices[v[2]] );
    lv0v1 = sqrt( ( Ov0 - Ov1 ).sq_length() );
    lv1v2 = sqrt( ( Ov1 - Ov2 ).sq_length() );
    lv2v0 = sqrt( ( Ov2 - Ov0 ).sq_length() );
    dv0 = min( lv0v1, lv2v0 );
    dv1 = min( lv1v2, lv0v1 );
    dv2 = min( lv2v0, lv1v2 );

    for(size_t i=0; i<3; i++){
        uint32_t ei = be[i];

        d = dv0;
        if(edges[ei].ep[1] == v[1])       d = dv1;
        else if(edges[ei].ep[1] == v[2])  d = dv2;
        d = min(d/3, epsilon);
        chamfer_edge_ep( ei, d, 1);

        d = dv0;
        if(edges[ei].ep[0] == v[1])       d = dv1;
        else if(edges[ei].ep[0] == v[2])  d = dv2;
        d = min(d/3, epsilon);
        chamfer_edge_ep( ei, d, 0);
    }
    
    remove_junk_edges_from_face(fi);

    for(size_t i=0; i < faces[fi].bounding_edges.size(); i++){
        if( edges[ faces[fi].bounding_edges[i] ].isAcute() ){
            uint32_t ei = faces[fi].bounding_edges[i];
            chamfer_acute_edge_from_inc_face( ei, fi, epsilon);
        }
    } 

}

// void PLCc::inputTriangleChamfering(uint32_t fi){
//     const CHAMface& f = faces[fi];

//     // Isolate f from incident faces:
//     // 1- remove fi from all edges bounding f
//     std::vector<std::pair<uint32_t,uint32_t>> rem_flat_edges;
//     for(const uint32_t& ei : f.bounding_edges) if(ei!=UINT32_MAX){ 
        
//         if(edges[ei].isFlat()){ 
//             rem_flat_edges.push_back(std::pair<uint32_t,uint32_t>(edges[ei].oep[0],edges[ei].oep[1]));
//             // DEBEUG
//             std::cout<<"WARNING detaching flat edge "<<ei<<" <"<<edges[ei].oep[0]<<","<<edges[ei].oep[1]<<"> while using fallback solution.\n";
//         }
//        else edges[ ei ].removeIncidentFace(fi);
//     }

//     // 2- substituite f with its related input triangle
//     uint32_t ti = f.triangle;
//     const uint32_t* in_tv = plc.triangle_vertices.data() + ti*3;
//     uint32_t v[] = { *in_tv, *(in_tv +1), *(in_tv +2) };
//     edges.push_back( CHAMedge(v[0], v[1], fi ) );
//     edges.push_back( CHAMedge(v[1], v[2], fi ) );
//     edges.push_back( CHAMedge(v[2], v[0], fi ) );
//     mark_edges.push_back(0);
//     mark_edges.push_back(0);
//     mark_edges.push_back(0);
    
//     faces[fi].bounding_edges.clear();
//     uint32_t l = (uint32_t) edges.size();
//     uint32_t be[] = {--l, --l, --l};
//     faces[fi].bounding_edges.assign( be, be+3 );
    
//     orient_initial_triface_bnd(faces[fi]);

//     // chamfer new face
//     edges[be[0]].type = CHAMedge_t::acute;
//     edges[be[1]].type = CHAMedge_t::acute;
//     edges[be[2]].type = CHAMedge_t::acute;

//     double lv0v1, lv1v2, lv2v0, dv0, dv1, dv2, d;
//     vector3d Ov0( vertices[v[0]] );
//     vector3d Ov1( vertices[v[1]] );
//     vector3d Ov2( vertices[v[2]] );
//     lv0v1 = sqrt( ( Ov0 - Ov1 ).sq_length() );
//     lv1v2 = sqrt( ( Ov1 - Ov2 ).sq_length() );
//     lv2v0 = sqrt( ( Ov2 - Ov0 ).sq_length() );
//     dv0 = min( lv0v1, lv2v0 );
//     dv1 = min( lv1v2, lv0v1 );
//     dv2 = min( lv2v0, lv1v2 );

//     for(size_t i=0; i<3; i++){
//         uint32_t ei = be[i];

//         d = dv0;
//         if(edges[ei].ep[1] == v[1])       d = dv1;
//         else if(edges[ei].ep[1] == v[2])  d = dv2;
//         d = min(d/3, epsilon);
//         chamfer_edge_ep( ei, d, 1);

//         d = dv0;
//         if(edges[ei].ep[0] == v[1])       d = dv1;
//         else if(edges[ei].ep[0] == v[2])  d = dv2;
//         d = min(d/3, epsilon);
//         chamfer_edge_ep( ei, d, 0);
//     }
    
//     remove_junk_edges_from_face(fi);

//     // identify flat edges
//     for(size_t i=0; i<3; i++){
//         uint32_t ei0 = edges[be[i]].oep[0];
//         uint32_t ei1 = edges[be[i]].oep[1];
//         for(const std::pair<uint32_t,uint32_t>& fep : rem_flat_edges){
//             if((fep.first == ei0 && fep.second == ei1) ||
//                (fep.first == ei1 && fep.second == ei0)   ){
//                     std::cout<<"recognized flat edge "<<be[i]<<" <"<<ei0<<","<<ei1<<"> while using fallback solution.\n";
//                     edges[be[i]].type == CHAMedge_t::flat;
//                 }
//         }
//     }

//     for(size_t i=0; i < faces[fi].bounding_edges.size(); i++){
//         if( edges[ faces[fi].bounding_edges[i] ].isAcute() && 
//             !edges[ faces[fi].bounding_edges[i] ].isFlat()    ){
//             uint32_t ei = faces[fi].bounding_edges[i];

//             // DEBUG
//             std::cout<<"edge "<<ei<<" sub-edge of <"<<edges[ei].oep[0]<<","<<edges[ei].oep[1]<<"> detached during fallback chamfering.\n";

//             chamfer_acute_edge_from_inc_face( ei, fi, epsilon);
//         }
//     } 

//     // DEBEUG
//     std::cout<<endl;

// }

// --------------------- //
// ORTHOGONAL CHAMFERING //
// --------------------- //

//
uint32_t PLCc::chamfering_face(uint32_t fi){

    #ifdef PLCC_VERBOSE_DEBUG
    std::cout<<"\nchamfering "; print_face_edges(fi);
    #endif

    // Part 1:
    // each couple of edges incident at an acute vertex will be
    // replaced by 3 new edges, i.e. the size of bounding_edges
    // have to be incresaed of the number of acute vertices on
    // face boundary.

    // Be shure that first and last bounding edges do not share an acute vertex
    uint32_t v0 = edges[ faces[fi].bounding_edges.front() ].commonVertex( edges[ faces[fi].bounding_edges.back() ] );
    if(is_acute_vrt( v0 )) faces[fi].make_first( faces[fi].bounding_edges.back() );

    size_t n_acute_vrts = 0;
    for(uint32_t& ei : faces[fi].bounding_edges) 
        if(is_acute_vrt(edges[ei].ep[0]) || is_acute_vrt(edges[ei].ep[1])) n_acute_vrts++;
    n_acute_vrts = n_acute_vrts / 2; // each vertex have been counted two times. 

    if(n_acute_vrts == 0) return 0;

    faces[fi].bounding_edges.resize( faces[fi].bounding_edges.size() + n_acute_vrts, EMPTY_PLACE);

    std::vector<uint32_t>& fbnd = faces[fi].bounding_edges;
    bool impossible = false;
    for(size_t i=0; i < fbnd.size()-1; i++) if(fbnd[i] != EMPTY_PLACE) {

        #ifdef PLCC_VERBOSE_DEBUG_LEV1
        std::cout<<"curr "; print_edge(fbnd[i]);
        #endif

        uint32_t ei = fbnd[i],  nei = fbnd[i+1];
        uint32_t vi = edges[ei].commonVertex( edges[nei] );

        if( is_acute_vrt(vi) ){
            // replace <v1,vi> + <vi,v4> with <v1,v2> + <v2,v3> + <v3,v4>

            #ifdef PLCC_VERBOSE_DEBUG
            std::cout<<"\nvertex "<<vi<<" between\n"; print_edge(ei); std::cout<<"and\n";
            print_edge(nei); std::cout<<"is acute.\n\n";
            #endif

            // 1) create two new vertices (v2 and v3) inside the face
            uint32_t v1, v4, ov1, ov4;
            v1 = edges[ei].oppositeVertex(vi);
            ov1 = edges[ei].oppositeOriginalVertex(vi);
            v4 = edges[nei].oppositeVertex(vi);
            ov4 = edges[nei].oppositeOriginalVertex(vi);
            uint32_t resp = new_vrts_in_inputTri(fi, vi, v1, ov1, v4, ov4);
            if(resp == INVALID_BPT){ impossible = true; break;}

            // 2) create 3 new edges
            uint32_t v2 = resp -1, v3 = resp;
            uint32_t num_e = edges.size()+3;
            edges.resize(num_e);
            edges[num_e-3].init_bridge_edge(v1, v2, vi, fi);
            edges[num_e-2].init_bridge_edge(v2, v3, vi, fi);
            edges[num_e-1].init_bridge_edge(v3, v4, vi, fi);
            mark_edges.resize(num_e, 0);

            #ifdef PLCC_DEBUG
            check_bridge(v1,v2,v3,v4);
            #endif
            #ifdef PLCC_VERBOSE_DEBUG
            std::cout<<"\nnew edges:\n"; 
            print_edge( (uint32_t) edges.size()-3 );
            print_edge( (uint32_t) edges.size()-2 );
            print_edge( (uint32_t) edges.size()-1 );
            #endif

            // 3) update boundary
            edges[ei].removeIncidentFace(fi);
            edges[nei].removeIncidentFace(fi);
            uint32_t pne = (uint32_t) edges.size()-1;
            fbnd[i] = pne-2; fbnd[i+1] = pne-1;
            for(size_t k=fbnd.size()-1; k>i+2; k--) fbnd[k] = fbnd[k-1];
            fbnd[i+2] = pne;

            #ifdef PLCC_VERBOSE_DEBUG
            std::cout<<"\nupdated face bnd "; print_face_edges(fi);
            #endif

            i += 2;
        }

    }

    if(impossible){ 

        // DEBUG
        const uint32_t* tv = plc.triangle_vertices.data() + faces[fi].triangle*3;
        std::cout<<"chamfering failed for face["<< fi <<"] (input tri #"
                 << faces[fi].triangle <<" <"<<*(tv)<<","<<*(tv+1)<<","<<*(tv+2)<<">) - fallback solution\n";
        
        if(verbose) std::cout<<"chamfering failed for face["<< fi <<"] -> fallback solution\n";

        // STUCTURAL-DEBUG
        for(const uint32_t& ei : faces[fi].bounding_edges){ 
            if(ei!=EMPTY_PLACE && edges[ei].isFlat()){
                std::cout<<"ERROR this face has the a sub edge of flat edge "<<ei<<" <"<<edges[ei].oep[0]<<","<<edges[ei].oep[1]<<"> on its boundary.\n";
               // ip_error("it is not possible to fallback chamfer - CLEVER SOLUTION NEED TO BE IMPLEMENTED\n");
               std::cout<<"[cham.cpp - chamfering_face()] ERROR Impossible to fallback chamfer\n"; exit(1);
            }
        }

        inputTriangleChamfering(fi); 

        return 1; 
    }

    // Part 2:
    // each acute edge will be "moved in to the face" by
    // replacing it and its two incident edges with a unique new edge.

    #ifdef PLCC_VERBOSE_DEBUG_LEV1
    std::cout<<"\nChamfer acute edges\n\n";
    #endif

    for(size_t i=0; i < fbnd.size(); i++ ) if( fbnd[i] != EMPTY_PLACE ){

        uint32_t ei = fbnd[i];

        #ifdef PLCC_VERBOSE_DEBUG_LEV1
        std::cout<<"curr "; print_edge(ei);
        #endif

        if( edges[ei].isAcute() ){

            assert( !edges[ei].isFlat() );

            #ifdef PLCC_VERBOSE_DEBUG
            std::cout<<"\nFound acute\n"; print_edge(ei);
            #endif
            
            const CHAMface& f = faces[fi];
            size_t ni = i, pi = i;
            f.advance_on_bnd(ni); while(fbnd[ni]==EMPTY_PLACE) f.advance_on_bnd(ni);
            f.reverse_on_bnd(pi); while(fbnd[pi]==EMPTY_PLACE) f.reverse_on_bnd(pi);
            uint32_t vn, vp;
            vn = edges[ fbnd[ni] ].notCommonVertex( edges[ei] );
            vp = edges[ fbnd[pi] ].notCommonVertex( edges[ei] );

            edges.push_back( CHAMedge( vn, vp, fi) );
            mark_edges.push_back(0);

            edges[ fbnd[ni] ].inc_face.clear();
            edges[ fbnd[pi] ].inc_face.clear();
            edges[ ei ].removeIncidentFace(fi);
            faces[fi].bounding_edges[i] = (uint32_t) edges.size()-1;
            faces[fi].bounding_edges[ni] = EMPTY_PLACE;
            faces[fi].bounding_edges[pi] = EMPTY_PLACE;

            #ifdef PLCC_VERBOSE_DEBUG
            std::cout<<"\nupdated face bnd "; print_face_edges(fi);
            #endif

            i = ni;
        }

    }

    CHAMface& f = faces[fi];
    f.bounding_edges.erase(
        std::remove_if( f.bounding_edges.begin(), 
                        f.bounding_edges.end(), 
                        [](uint32_t x){ return x == EMPTY_PLACE; } ), 
        f.bounding_edges.end());

    #ifdef PLCC_DEBUG
    std::cout<<"after removing..\n"; print_face_edges(fi);
    #endif

    return 0;
}

// --------------- //
// CHAMFERING MAIN //
// --------------- //

void PLCc::chamfering(){

    chamfering_vrts();
    if(verbose) std::cout<<"[chamfering()] vertex chamfering COMPLETED\n";

    #ifdef PLCC_DEBUG
    checkEdges_beforeJunkDeletion();
    #endif

    uint32_t n_fbs = 0;
    for(uint32_t fi=0; fi<faces.size(); fi++){ n_fbs += chamfering_face(fi); }

    // Delete junk edges (and all isoleted edges) from edges vector
    for(CHAMedge& e : edges) if( e.isJunk() ){ e.inc_face.clear(); }
    cleanUp_edges();

    if(n_fbs) std::cout<<n_fbs<<" (on "<<faces.size()<<") faces have be fall-back chamfered.\n";
    if(verbose) std::cout<<"[chamfering()] edge chamfering COMPLETED\n";
}

// ---------------------------------------------- //
// CHAMFERED PLC SIMPLIFICATION (Post-Processing) //
// ---------------------------------------------- //

class bridge{
    public:
    uint32_t el, ec, er;
    double shortest_e_sq;

    uint32_t s, sl, sr;

    bridge() : el(UINT32_MAX), ec(UINT32_MAX), er(UINT32_MAX), shortest_e_sq(DBL_MAX),
               s(UINT32_MAX), sl(UINT32_MAX), sr(UINT32_MAX) {}

    inline bool is_3bridge() const { return el!=UINT32_MAX && er!=UINT32_MAX; }
    inline bool is_2bridge() const { return (el==UINT32_MAX && er!= UINT32_MAX) || (el!=UINT32_MAX && er== UINT32_MAX); }

    static inline bool bridgeSortFunction(const bridge& b1, const bridge& b2) { return (b1.shortest_e_sq < b2.shortest_e_sq); }
};

inline std::ostream& operator<<(std::ostream& os, const bridge& b) {
    os <<"<";
    if(b.el!=UINT32_MAX) os << (b.el); else os << "empty";
     os << ", " << b.ec << ", " ;
    if(b.er!=UINT32_MAX) os << (b.er); else os << "empty";
    os << "> - short_len = " << b.shortest_e_sq;
    if( b.s != UINT32_MAX ) os << " => <" << b.s << ">";
    else if( b.sl != UINT32_MAX ) os <<  " => <" << b.sl << ", "<< b.sr <<">";
    return os;
}

// NOTE: Assumes that edges[ei] has a unique incident face,
//       AND that conn_edge has exactly two incident faces.
// In practise ei, conn_edge and cons_edge forms meets in a unique vertex.
uint32_t PLCc::get_cons_edge_on_adj_face(uint32_t ei){
    uint32_t conn_edge, cons_edge;
    uint32_t fi = edges[ei].inc_face[0];
    faces[fi].make_first(ei); 
    const std::vector<uint32_t>& fbnd = faces[fi].bounding_edges;
    conn_edge = fbnd[1];
    if(edges[conn_edge].inc_face.size() == 1) conn_edge = fbnd.back();
    if(edges[conn_edge].inc_face.size() != 2) return UINT32_MAX;
    uint32_t fj = edges[conn_edge].inc_face[0];
    if(fj == fi) fj = edges[conn_edge].inc_face[1];
    faces[fj].make_first(conn_edge);
    cons_edge = faces[fj].bounding_edges[1];
    if(!edges[cons_edge].has_commonVertex(edges[ei])) cons_edge = faces[fj].bounding_edges.back();
    assert(edges[cons_edge].has_commonVertex(edges[ei]));
    return cons_edge;
}

void PLCc::chamfered_plc_simplification(){

    // Collect bridge-edges
    std::vector<bridge> bridges;
    uint32_t e0, e1, ei;
    for(CHAMface& f : faces) {
        for(uint32_t fb : f.bounding_edges) if(edges[fb].loc_face_bridge_id == NO_BRIDGE) f.make_last(fb);
        const std::vector<uint32_t>& fbnd = f.bounding_edges;

        for(size_t i=0; i < (uint32_t) fbnd.size(); i++) if(edges[ fbnd[i] ].loc_face_bridge_id != NO_BRIDGE) {
            bridges.push_back( bridge() );
            bridge& b = bridges.back();
            do{ 
                ei = fbnd[i];
                const CHAMedge& e = edges[ei];
                e0 = e.ep[0]; e1 = e.ep[1];
                if( vertices[ e0 ]->isLNC() && vertices[ e1 ]->isBPT() ) b.el = ei;
                else if( vertices[ e0 ]->isBPT() && vertices[ e1 ]->isBPT() ) b.ec = ei;
                else if( vertices[ e0 ]->isBPT() && vertices[ e1 ]->isLNC() ) b.er = ei;
                else {
                    // ip_error("[chamfered_plc_simplification] invalid bridge-edge\n");
                    std::cout<<"[cham.cpp - chamfered_plc_simplification()] ERROR Invalid bridge-edge\n"; exit(1);
                }
                b.shortest_e_sq = min( eEdgeSqLen(ei), b.shortest_e_sq );
                i++;
            }while(edges[ fbnd[i] ].loc_face_bridge_id == edges[ fbnd[i-1] ].loc_face_bridge_id);
            if(edges[ fbnd[i] ].loc_face_bridge_id!=NO_BRIDGE) i--;
        }
    }

    if(verbose){ 
        uint32_t nbrs = 0;  
        for(const CHAMedge& e : edges){ if(e.loc_face_bridge_id != NO_BRIDGE) nbrs++; }
        if(nbrs > 0) std::cout<<"There are "<< nbrs << " bridge-edges and "<< bridges.size()<< " bridges.\n";
    }

    #ifdef PLCC_DEBUG_VERBOSE_LEV1
    print_all_bridge_edges();
    #endif 
    #ifdef PLCC_DEBUG
    check_bridge_edges();
    #endif

    // Order bridges by increasing shortest length
    std::sort( bridges.begin(), bridges.end(), bridge::bridgeSortFunction );

    // DEBUG
    // std::cout<<"\n";
    // for(const bridge& b : bridges) std::cout<<b<<"\n";
    // std::cout<<"\n";

    // Simplify bridges if acute angles are not introducted
    for(size_t i=0; i<bridges.size(); i++){
        bridge& b = bridges[i];

        // build a unique edge alternative
        if( b.is_3bridge() ){
            uint32_t cons_left_edge = get_cons_edge_on_adj_face(b.el);
            uint32_t cons_right_edge = get_cons_edge_on_adj_face(b.er);

            // DEBUG
            // std::cout<<"\nleft_ext: "; print_edge(cons_left_edge);
            // std::cout<<"left_bridge: "; print_edge(b.el);  
            // std::cout<<"cent_bridge: "; print_edge(b.ec);
            // std::cout<<"right_bridge: "; print_edge(b.er);
            // std::cout<<"right_ext: "; print_edge(cons_right_edge);

            // DEBUG
            // if( cons_left_edge==UINT32_MAX || cons_right_edge==UINT32_MAX ){

            //     if( cons_left_edge==UINT32_MAX )
            //         std::cout<<"ERROR: cons_left_edge NOT FOUND!\n";
            //      if( cons_right_edge==UINT32_MAX )  
            //         std::cout<<"ERROR: cons_right_edge NOT FOUND!\n";

            //     std::cout<<"left_bridge: "; print_edge(b.el);  
            //     std::cout<<"cent_bridge: "; print_edge(b.ec);
            //     std::cout<<"right_bridge: "; print_edge(b.er);

            //     assert(edges[ b.ec ].inc_face.size() == 1);
            //     const CHAMface& f = faces[ edges[ b.ec ].inc_face[0] ];
            //     for(size_t bei = 0; bei < f.bounding_edges.size(); bei++){
            //         print_edge_and_inc_face(f.bounding_edges[bei]);
            //     }
            // }

            bool simpL = true, simpR = true;
            uint32_t v0, v1, v2, v3;
            v1 = edges[b.el].notCommonVertex(edges[b.ec]);
            v2 = edges[b.er].notCommonVertex(edges[b.ec]);
            if(cons_left_edge==UINT32_MAX) simpL = false;
            else{
                v0 = edges[cons_left_edge].notCommonVertex(edges[b.el]);
                simpL = !isAcuteAngle(vertices[v0], vertices[v1], vertices[v2]);
            }
            if(cons_right_edge==UINT32_MAX) simpR = false;
            else{
                v3 = edges[cons_right_edge].notCommonVertex(edges[b.er]);
                simpR = !isAcuteAngle(vertices[v1], vertices[v2], vertices[v3]);
            }

            // DEBUG
            // uint32_t vc0 = edges[b.el].commonVertex(edges[b.ec]);
            // uint32_t vc1 = edges[b.er].commonVertex(edges[b.ec]);
            // bool simpLC = !isAcuteAngle(vertices[v1], vertices[vc0], vertices[vc1]);
            // bool simpRC = !isAcuteAngle(vertices[vc0], vertices[vc1], vertices[v2]);
            // std::cout<<"CL: "<<simpLC<<" , CR: "<<simpRC<<"\n";
            // std::cout<<"L: "<<simpL<<" , R: "<<simpR<<"\n\n";

            if(simpL && simpR){
                uint32_t fi = edges[b.ec].inc_face[0];
                edges.push_back( CHAMedge(v1,v2,fi) ); mark_edges.push_back(0);
                size_t new_fbnd_size = faces[fi].bounding_edges.size()-2 ;
                faces[ fi ].make_last(b.er);
                if(faces[fi].bounding_edges[ new_fbnd_size ] != b.ec) faces[ fi ].make_last(b.el);
                faces[fi].bounding_edges[ new_fbnd_size-1 ] = (uint32_t)edges.size()-1; 
                faces[fi].bounding_edges.resize(new_fbnd_size);
                edges[b.el].isolate();
                edges[b.ec].isolate();
                edges[b.er].isolate();
            }
            else if(simpL || simpR){
                uint32_t fi = edges[b.ec].inc_face[0];
                uint32_t be_rem, new_ep, ext_ep; 
                if(simpR){ be_rem = b.er; new_ep = v2; ext_ep = v3;}
                else{ be_rem = b.el, new_ep = v1, ext_ep = v0; }
                uint32_t c = edges[b.ec].notCommonVertex(edges[be_rem]);
                uint32_t d = edges[b.ec].commonVertex(edges[be_rem]);

                // Further simplification: same lenght (TO IMPLEMENT)
                // double curr_btp_u = vertices[c]->toBPT().U();
                // double curr_btp_v = vertices[c]->toBPT().V();
                // double optim_bpt_u = (curr_btp_u + vertices[d]->toBPT().U()) * 0.5;
                // double optim_bpt_v = (curr_btp_v + vertices[d]->toBPT().V()) * 0.5;
                // const implicitPoint3D_BPT* candidate = new implicitPoint3D_BPT( vertices[c]->toBPT().P(), vertices[c]->toBPT().Q(), vertices[c]->toBPT().R(), optim_bpt_u, optim_bpt_v );

                if(!isAcuteAngle(vertices[v1],vertices[c],vertices[v2]) &&
                   !isAcuteAngle(vertices[ext_ep],vertices[new_ep],vertices[c]) ){

                    edges.push_back( CHAMedge(new_ep,c,fi) ); mark_edges.push_back(0);
                    faces[ fi ].make_last(be_rem);
                    size_t new_fbnd_size = faces[fi].bounding_edges.size()-1;
                    if( faces[fi].bounding_edges[ new_fbnd_size-1 ] != b.ec){ 
                        faces[ fi ].make_last(b.ec);
                        assert(faces[fi].bounding_edges[ new_fbnd_size-1 ] == be_rem);
                    }
                    faces[fi].bounding_edges[ new_fbnd_size-1 ] = (uint32_t)edges.size()-1; 
                    faces[fi].bounding_edges.resize(new_fbnd_size);
                    edges[b.ec].isolate();
                    edges[be_rem].isolate();
                }
            }
        } 
        else if(b.is_2bridge()){
            uint32_t cons_edge;
            if(b.er==UINT32_MAX) cons_edge = get_cons_edge_on_adj_face(b.el);
            else cons_edge = get_cons_edge_on_adj_face(b.er);

            // // DEBUG
            // std::cout<<"\next: "; print_edge(cons_edge);
            // std::cout<<"side_bridge: "; if(b.er==UINT32_MAX) print_edge(b.el); else print_edge(b.er);
            // std::cout<<"cent_bridge: "; print_edge(b.ec);

            if(cons_edge!=UINT32_MAX){

                uint32_t be = b.er;
                if(be==UINT32_MAX) be = b.el;
                uint32_t v0 = edges[cons_edge].notCommonVertex(edges[be]);
                uint32_t v1 = edges[cons_edge].commonVertex(edges[be]);
                uint32_t v2 = edges[b.ec].notCommonVertex(edges[be]);

                if(!isAcuteAngle(vertices[v0], vertices[v1], vertices[v2])){
                    uint32_t fi = edges[b.ec].inc_face[0];
                    edges.push_back( CHAMedge(v1,v2,fi) ); mark_edges.push_back(0);
                    faces[ fi ].make_last(be);
                    size_t new_fbnd_size = faces[fi].bounding_edges.size()-1 ;
                    if( faces[fi].bounding_edges[ new_fbnd_size-1 ] != b.ec){ 
                        faces[ fi ].make_last(b.ec);
                        assert(faces[fi].bounding_edges[ new_fbnd_size-1 ] == be);
                    }
                    faces[fi].bounding_edges[ new_fbnd_size-1 ] = (uint32_t)edges.size()-1; 
                    faces[fi].bounding_edges.resize(new_fbnd_size);
                    edges[b.ec].isolate();
                    edges[be].isolate();
                }

                
            }

        }

    }

    // Remove superfluous edges
    cleanUp_edges();
}


// ------------------- //
// FINAL TRIANGULATION //
// ------------------- //

//
bool isConvexVrt(const pointType *p, const pointType *vrt, const pointType *q, int* ref_nv_sign)
{
    // DEBUG
    // std::cout<<"tri norm vect signs ("<<ref_nv_sign[0]<<","<<
    //                                     ref_nv_sign[1]<<","<<
    //                                     ref_nv_sign[2]<<")\n";
    //  std::cout<<"loc norm vect signs ("<<pointType::orient2Dyz(*p,*vrt,*q)<<","<<
    //                                     pointType::orient2Dzx(*p,*vrt,*q)<<","<<
    //                                     pointType::orient2Dxy(*p,*vrt,*q)<<")\n";
    
    if(pointType::orient2Dyz(*p,*vrt,*q) != ref_nv_sign[0]) return false;
    if(pointType::orient2Dzx(*p,*vrt,*q) != ref_nv_sign[1]) return false;
    if(pointType::orient2Dxy(*p,*vrt,*q) != ref_nv_sign[2]) return false;
    return true;
}

//
bool pointInTriangle(const pointType* pt, const pointType* v0, const pointType* v1, const pointType* v2){
    // Returns true if the coplanar point 'p' is either in the inner area of
    // 't' or on its border. Undetermined if p and t are not coplanar.
	if (pointType::pointInSegment(*pt, *v0, *v1)) return true;
	if (pointType::pointInSegment(*pt, *v1, *v2)) return true;
	if (pointType::pointInSegment(*pt, *v2, *v0)) return true;
	return pointType::pointInInnerTriangle(*pt, *v0, *v1, *v2);
}

// HEAR CLIPPING TRIANGULATION ALGORITHM
// Triangulation of a 2D polygon.
// Assumptions: the polygon is simple (no holes, closed, no self-intersections).
void PLCc::hear_clipping(uint32_t fi, std::vector<uint32_t>& out_tri_fv_list) const {

    const uint32_t* in_tri = plc.triangle_vertices.data() + (faces[fi].triangle*3);
    const pointType* ov0 = vertices[*in_tri];
    const pointType* ov1 = vertices[*(in_tri+1)];
    const pointType* ov2 = vertices[*(in_tri+2)];
    int signs_normVect_tri[] = {pointType::orient2Dyz(*ov0,*ov1,*ov2),
                                pointType::orient2Dzx(*ov0,*ov1,*ov2),
                                pointType::orient2Dxy(*ov0,*ov1,*ov2) };
            
    // Find first ear vertex. Create triangle. Remove ear vertex from list.
    // Do this while number of vertices (fv.size()) > 2.

    std::vector<uint32_t> fv; get_face_vertices(faces[fi], fv);
    size_t jl, jc = 0, jr;
    size_t fv_last = fv.size()-1;
    while(fv.size() > 2){

        if(jc == 0){ jl = fv_last; jr = 1; }
        else{ 
            if(jc == fv_last){ jl = jc-1; jr = 0; }
            else{ jl = jc-1; jr = jc+1; }
        }

        // DEBUG
        // std::cout<<"\njl = "<<jl<<", jc = "<<jc<<", jr = "<<jr<<"\n";
        // std::cout<<" fv["<<jl<<"] = "<<fv[jl]<<
        //            ", fv["<<jc<<"] = "<<fv[jc]<<
        //            ", fv["<<jr<<"] = "<<fv[jr]<<"\n";

        const pointType* vl = vertices[ fv[jl] ];
        const pointType* vc = vertices[ fv[jc] ];
        const pointType* vr = vertices[ fv[jr] ];

        bool is_convex = isConvexVrt(vl,vc,vr,signs_normVect_tri);

        // DEBUG
        // std::cout<<fv[jc]<<" is ";
        // if(!is_convex) std::cout<<"NOT ";
        // std::cout<<"convex vertex.\n";

        bool is_ear = true;

        if(is_convex){

            size_t jtest = 0;
            if(jr != fv_last) jtest = jr+1;

            // DEBUG
            // std::cout<<"jtest = "<<jtest<<"; fv["<<jtest<<"] = "<<fv[jtest]<<"\n";

            while( jtest != jl && is_ear ){
                const pointType* vt = vertices[ fv[jtest] ];
                is_ear = !pointInTriangle(vt,vl,vc,vr);

                // DEBUG
                // std::cout<<fv[jtest]<<" is ";
                // if(is_ear) std::cout<<"NOT ";
                // std::cout<<"inside <"<<fv[jl]<<","<<fv[jc]<<","<<fv[jr]<<">\n";

                if(jtest==fv_last) jtest = 0; else jtest++;

                // DEBUG
                // if(is_ear) std::cout<<"next jtest = "<<jtest<<"; fv["<<jtest<<"] = "<<fv[jtest]<<"\n";
            }

        }
        else{ is_ear = false; }

        if(is_ear){

            // DEBUG
            // std::cout<<"detraching triangle <"<<fv[jl]<<","<<fv[jc]<<","<<fv[jr]<<">\n";

            out_tri_fv_list.insert(out_tri_fv_list.end(), {fv[jl],fv[jc],fv[jr]} );
            fv.erase(std::remove(fv.begin(),fv.end(),fv[jc]));
            fv_last = fv.size()-1;
            jc = 0;
        }
        else{
            if(jc==fv_last) jc = 0; else jc++;
        }

    }
}

// 
void PLCc::get_triangles(std::vector<uint32_t>& tri_fv) const {

    for(size_t fi = 0; fi < faces.size(); fi++) {

        // std::cout<<"\ntriangualting "; print_face_edges(fi); // DEBUG

        if(faces[fi].bounding_edges.size() == 3){
            std::vector<uint32_t> fv;
            get_face_vertices(faces[fi], fv);
            tri_fv.insert(tri_fv.end(), {fv[0],fv[1],fv[2]} );
        }
        else{

            // std::cout<<"HEAR CLIPPING\n"; // DEBUG

            hear_clipping((uint32_t)fi, tri_fv);
        }

    }
}

// Complmentar
void PLCc::get_complementar_tri(const std::vector<uint32_t>& out_tri, std::vector<uint32_t>& compl_tri){
    // make half-edges
    std::vector< std::vector<uint32_t> > he; // ep0, ep1, tri, occurences for each edge
    uint32_t v0,v1,v2,e0,e1;
    for(uint32_t i=0; i<out_tri.size()/3; i++){
        v0 = out_tri[3*i  ];
        v1 = out_tri[3*i+1];
        v2 = out_tri[3*i+2];
        e0 = v0; e1 = v1; if(v0 > v1) {e0 = v1; e1 = v0;} he.push_back({e0,e1,i,1});
        e0 = v1; e1 = v2; if(v1 > v2) {e0 = v2; e1 = v1;} he.push_back({e0,e1,i,1});
        e0 = v2; e1 = v0; if(v2 > v0) {e0 = v0; e1 = v2;} he.push_back({e0,e1,i,1});
    }
    sort(he.begin(), he.end(),
        [](const std::vector<uint32_t> &a, const std::vector<uint32_t> &b){ return (a[0] < b[0] || (a[0]==b[0] && a[1]<b[1])); } );

    for(size_t i=0; i<he.size()-1; i++){
        if(he[i][0] == he[i+1][0] && he[i][1] == he[i+1][1]){
            he[i+1][3] += he[i][3];
            he[i][3] = 0;
        }
    }

    he.erase(std::remove_if(he.begin(), he.end(), [](const std::vector<uint32_t> &a){ return (a[3] != 1); }), he.end());
    
    uint32_t r0, r1;
    for(size_t i=0; i<he.size(); i++){
        e0 = he[i][0]; e1 = he[i][1]; 

        if(vertices[e0]->isExplicit3D() || vertices[e1]->isExplicit3D() ) continue;

        r0 = ref_exp3D_vrt[e0]; r1 = ref_exp3D_vrt[e1];
        if( r0 == r1 ) {
            compl_tri.insert( compl_tri.end(), {e0, e1, r0} ); // orientation maybe have to be changed
        }
        else{
            if(e0 > e1) std::swap(e0, e1);
            compl_tri.insert( compl_tri.end(), {e0, e1, r0} );
            compl_tri.insert( compl_tri.end(), {r0, r1, e1} );
            // orientation maybe have to be changed
        }
    }

    for(size_t i=0; i<compl_tri.size()/3; i++){
        uint32_t a = compl_tri[3*i];
        uint32_t b = compl_tri[3*i+1];
        uint32_t c = compl_tri[3*i+2];
        assert( !vCollinear(a,b,c) );
    }
}

// WORK IN PROGRESS
void PLCc::triangulate_chamfered_plc(){
    std::cout<<"[triangulate_chamfered_plc] THIS FUNCTION IS NOT COMPLETE.";
    assert( checkup() );
}


// To export/print a triangulated output without modifing the data structure

void PLCc::get_trivrts_1av(const std::vector<uint32_t>& fv, std::vector<uint32_t>& tri_fv) const {

    const uint32_t nv = 6;
    assert(fv.size() == nv);

    uint32_t EXP1, LNC1, BPT1, BPT2, LNC2, EXP2;
    uint32_t i_ref;
    for (size_t i = 0; i < nv; i++) {
        uint32_t i1 = (i + 1) % nv;
        if (vertices[fv[i]]->isExplicit3D() &&
            vertices[fv[i1]]->isExplicit3D()) {
            i_ref = i1;
            break;
        }
    }

    EXP1 = fv[(i_ref + 0) % nv]; assert(vertices[EXP1]->isExplicit3D());
    LNC1 = fv[(i_ref + 1) % nv]; assert(vertices[LNC1]->isLNC());
    BPT1 = fv[(i_ref + 2) % nv]; assert(vertices[BPT1]->isBPT());
    BPT2 = fv[(i_ref + 3) % nv]; assert(vertices[BPT2]->isBPT());
    LNC2 = fv[(i_ref + 4) % nv]; assert(vertices[LNC2]->isLNC());
    EXP2 = fv[(i_ref + 5) % nv]; assert(vertices[EXP2]->isExplicit3D());

    // bool variant = false;
    // if( vInnerSegmentCross(EXP1, BPT2, LNC1, BPT1) ){ 
    //     assert( !vInnerSegmentCross(EXP2, BPT1, LNC2, BPT2) );
    //     variant = true;
    // }

    tri_fv.insert(tri_fv.end(), { EXP1, LNC1, BPT1 });
    tri_fv.insert(tri_fv.end(), { EXP1, BPT1, BPT2 });
    tri_fv.insert(tri_fv.end(), { EXP1, BPT2, EXP2 });
    // if(!variant){
    //     tri_fv.insert( tri_fv.end(), {EXP1, BPT1, BPT2} );
    //     tri_fv.insert( tri_fv.end(), {EXP1, BPT2, EXP2} );
    // }
    // else{
    //     tri_fv.insert( tri_fv.end(), {EXP1, BPT1, EXP2} );
    //     tri_fv.insert( tri_fv.end(), {EXP2, BPT1, BPT2} );
    // }
    tri_fv.insert(tri_fv.end(), { EXP2, BPT2, LNC2 });
}

void PLCc::get_trivrts_2av(const std::vector<uint32_t>& fv, std::vector<uint32_t>& tri_fv) const {

    const uint32_t nv = 9;
    assert(fv.size() == nv);

    uint32_t EXP, LNC1, BPT1, BPT2, LNC2, LNC3, BPT3, BPT4, LNC4;
    uint32_t i_ref;
    for (uint32_t i = 0; i < nv; i++) {
        if (vertices[fv[i]]->isExplicit3D()) {
            i_ref = i;
            break;
        }
    }

    EXP = fv[(i_ref + 0) % nv]; assert(vertices[EXP]->isExplicit3D());
    LNC1 = fv[(i_ref + 1) % nv]; assert(vertices[LNC1]->isLNC());
    BPT1 = fv[(i_ref + 2) % nv]; assert(vertices[BPT1]->isBPT());
    BPT2 = fv[(i_ref + 3) % nv]; assert(vertices[BPT2]->isBPT());
    LNC2 = fv[(i_ref + 4) % nv]; assert(vertices[LNC2]->isLNC());
    LNC3 = fv[(i_ref + 5) % nv]; assert(vertices[LNC3]->isLNC());
    BPT3 = fv[(i_ref + 6) % nv]; assert(vertices[BPT3]->isBPT());
    BPT4 = fv[(i_ref + 7) % nv]; assert(vertices[BPT4]->isBPT());
    LNC4 = fv[(i_ref + 8) % nv]; assert(vertices[LNC4]->isLNC());

    tri_fv.insert(tri_fv.end(), { EXP, LNC1, BPT1 });
    tri_fv.insert(tri_fv.end(), { EXP, BPT1, BPT2 });
    tri_fv.insert(tri_fv.end(), { EXP, BPT2, BPT3 });
    tri_fv.insert(tri_fv.end(), { EXP, BPT3, BPT4 });
    tri_fv.insert(tri_fv.end(), { EXP, BPT4, LNC4 });
    tri_fv.insert(tri_fv.end(), { BPT2, LNC2, BPT3 });
    tri_fv.insert(tri_fv.end(), { BPT3, LNC2, LNC3 });
}

void PLCc::get_trivrts_2av1ae(const std::vector<uint32_t>& fv, std::vector<uint32_t>& tri_fv) const {

    const uint32_t nv = 7;
    assert(fv.size() == nv);

    uint32_t EXP, LNC1, BPT1, BPT2, BPT3, BPT4, LNC2;
    uint32_t i_ref;
    for (uint32_t i = 0; i < nv; i++) {
        if (vertices[fv[i]]->isExplicit3D()) {
            i_ref = i;
            break;
        }
    }

    EXP = fv[(i_ref + 0) % nv];  assert(vertices[EXP]->isExplicit3D());
    LNC1 = fv[(i_ref + 1) % nv];  assert(vertices[LNC1]->isLNC());
    BPT1 = fv[(i_ref + 2) % nv];  assert(vertices[BPT1]->isBPT());
    BPT2 = fv[(i_ref + 3) % nv];  assert(vertices[BPT2]->isBPT());
    BPT3 = fv[(i_ref + 4) % nv];  assert(vertices[BPT3]->isBPT());
    BPT4 = fv[(i_ref + 5) % nv];  assert(vertices[BPT4]->isBPT());
    LNC2 = fv[(i_ref + 6) % nv];  assert(vertices[LNC2]->isLNC());

    tri_fv.insert(tri_fv.end(), { EXP, LNC1, BPT1 });
    tri_fv.insert(tri_fv.end(), { EXP, BPT1, BPT2 });
    tri_fv.insert(tri_fv.end(), { EXP, BPT2, BPT3 });
    tri_fv.insert(tri_fv.end(), { EXP, BPT3, BPT4 });
    tri_fv.insert(tri_fv.end(), { EXP, BPT4, LNC2 });
}

void PLCc::get_trivrts_3av(const std::vector<uint32_t>& fv, std::vector<uint32_t>& tri_fv) const {

    const uint32_t nv = 12;
    assert(fv.size() == nv);

    uint32_t LNC1, BPT1, BPT2, LNC2, LNC3, BPT3, BPT4, LNC4, LNC5, BPT5, BPT6, LNC6;
    uint32_t i_ref;
    for (uint32_t i = 0; i < nv; i++) {
        uint32_t i1 = (i + 1) % nv;
        if (vertices[fv[i]]->isLNC() &&
            vertices[fv[i1]]->isBPT()) {
            i_ref = i;
            break;
        }
    }

    LNC1 = fv[(i_ref + 0) % nv];  assert(vertices[LNC1]->isLNC());
    BPT1 = fv[(i_ref + 1) % nv];  assert(vertices[BPT1]->isBPT());
    BPT2 = fv[(i_ref + 2) % nv];  assert(vertices[BPT2]->isBPT());
    LNC2 = fv[(i_ref + 3) % nv];  assert(vertices[LNC2]->isLNC());
    LNC3 = fv[(i_ref + 4) % nv];  assert(vertices[LNC3]->isLNC());
    BPT3 = fv[(i_ref + 5) % nv];  assert(vertices[BPT3]->isBPT());
    BPT4 = fv[(i_ref + 6) % nv];  assert(vertices[BPT4]->isBPT());
    LNC4 = fv[(i_ref + 7) % nv];  assert(vertices[LNC4]->isLNC());
    LNC5 = fv[(i_ref + 8) % nv];  assert(vertices[LNC5]->isLNC());
    BPT5 = fv[(i_ref + 9) % nv];  assert(vertices[BPT5]->isBPT());
    BPT6 = fv[(i_ref + 10) % nv]; assert(vertices[BPT6]->isBPT());
    LNC6 = fv[(i_ref + 11) % nv]; assert(vertices[LNC6]->isLNC());

    tri_fv.insert(tri_fv.end(), { BPT1, BPT3, BPT5 });
    tri_fv.insert(tri_fv.end(), { BPT1, BPT2, BPT3 });
    tri_fv.insert(tri_fv.end(), { BPT3, BPT4, BPT5 });
    tri_fv.insert(tri_fv.end(), { BPT5, BPT6, BPT1 });
    tri_fv.insert(tri_fv.end(), { BPT6, LNC1, BPT1 });
    tri_fv.insert(tri_fv.end(), { BPT6, LNC6, LNC1 });
    tri_fv.insert(tri_fv.end(), { BPT2, LNC3, BPT3 });
    tri_fv.insert(tri_fv.end(), { BPT2, LNC2, LNC3 });
    tri_fv.insert(tri_fv.end(), { BPT4, LNC4, BPT5 });
    tri_fv.insert(tri_fv.end(), { BPT5, LNC4, LNC5 });
}

void PLCc::get_trivrts_3av1ae(const std::vector<uint32_t>& fv, std::vector<uint32_t>& tri_fv) const {

    const uint32_t nv = 10;
    assert(fv.size() == nv);

    uint32_t LNC1, BPT1, BPT2, LNC2, LNC3, BPT3, BPT4, BPT5, BPT6, LNC4;
    uint32_t i_ref;
    for (uint32_t i = 0; i < nv; i++) {
        uint32_t i1 = (i + 1) % nv;
        uint32_t i2 = (i + 2) % nv;
        uint32_t i3 = (i + 3) % nv;
        uint32_t i4 = (i + 4) % nv;
        uint32_t i5 = (i + 5) % nv;
        if (vertices[fv[i]]->isBPT() &&
            vertices[fv[i1]]->isBPT() &&
            vertices[fv[i2]]->isBPT() &&
            vertices[fv[i3]]->isBPT() &&
            vertices[fv[i4]]->isLNC()) {
            i_ref = i5;
            break;
        }
    }

    LNC1 = fv[(i_ref + 0) % nv]; assert(vertices[LNC1]->isLNC());
    BPT1 = fv[(i_ref + 1) % nv]; assert(vertices[BPT1]->isBPT());
    BPT2 = fv[(i_ref + 2) % nv]; assert(vertices[BPT2]->isBPT());
    LNC2 = fv[(i_ref + 3) % nv]; assert(vertices[LNC2]->isLNC());
    LNC3 = fv[(i_ref + 4) % nv]; assert(vertices[LNC3]->isLNC());
    BPT3 = fv[(i_ref + 5) % nv]; assert(vertices[BPT3]->isBPT());
    BPT4 = fv[(i_ref + 6) % nv]; assert(vertices[BPT4]->isBPT());
    BPT5 = fv[(i_ref + 7) % nv]; assert(vertices[BPT5]->isBPT());
    BPT6 = fv[(i_ref + 8) % nv]; assert(vertices[BPT6]->isBPT());
    LNC4 = fv[(i_ref + 9) % nv]; assert(vertices[LNC4]->isLNC());

    tri_fv.insert(tri_fv.end(), { BPT1, BPT3, BPT5 });
    tri_fv.insert(tri_fv.end(), { BPT1, BPT2, BPT3 });
    tri_fv.insert(tri_fv.end(), { BPT3, BPT4, BPT5 });
    tri_fv.insert(tri_fv.end(), { BPT5, BPT6, BPT1 });
    tri_fv.insert(tri_fv.end(), { BPT6, LNC1, BPT1 });
    tri_fv.insert(tri_fv.end(), { BPT6, LNC4, LNC1 });
    tri_fv.insert(tri_fv.end(), { BPT2, LNC3, BPT3 });
    tri_fv.insert(tri_fv.end(), { BPT2, LNC2, LNC3 });
}

void PLCc::get_trivrts_3av2ae(const std::vector<uint32_t>& fv, std::vector<uint32_t>& tri_fv) const {

    const uint32_t nv = 8;
    assert(fv.size() == nv);

    uint32_t LNC1, BPT1, BPT2, BPT3, BPT4, BPT5, BPT6, LNC2;
    uint32_t i_ref;
    for (uint32_t i = 0; i < nv; i++) {
        uint32_t i1 = (i + 1) % nv;
        if (vertices[fv[i]]->isLNC() &&
            vertices[fv[i1]]->isLNC()) {
            i_ref = i1;
            break;
        }
    }

    LNC1 = fv[(i_ref + 0) % nv]; assert(vertices[LNC1]->isLNC());
    BPT1 = fv[(i_ref + 1) % nv]; assert(vertices[BPT1]->isBPT());
    BPT2 = fv[(i_ref + 2) % nv]; assert(vertices[BPT2]->isBPT());
    BPT3 = fv[(i_ref + 3) % nv]; assert(vertices[BPT3]->isBPT());
    BPT4 = fv[(i_ref + 4) % nv]; assert(vertices[BPT4]->isBPT());
    BPT5 = fv[(i_ref + 5) % nv]; assert(vertices[BPT5]->isBPT());
    BPT6 = fv[(i_ref + 6) % nv]; assert(vertices[BPT6]->isBPT());
    LNC2 = fv[(i_ref + 7) % nv]; assert(vertices[LNC2]->isLNC());

    tri_fv.insert(tri_fv.end(), { BPT1, BPT3, BPT5 });
    tri_fv.insert(tri_fv.end(), { BPT1, BPT2, BPT3 });
    tri_fv.insert(tri_fv.end(), { BPT3, BPT4, BPT5 });
    tri_fv.insert(tri_fv.end(), { BPT5, BPT6, BPT1 });
    tri_fv.insert(tri_fv.end(), { BPT6, LNC1, BPT1 });
    tri_fv.insert(tri_fv.end(), { BPT6, LNC2, LNC1 });
}

void PLCc::get_trivrts_3av3ae(const std::vector<uint32_t>& fv, std::vector<uint32_t>& tri_fv) const {

    assert(fv.size() == 6);

    for (uint32_t v : fv) assert(vertices[v]->isBPT());

    tri_fv.insert(tri_fv.end(), { fv[0], fv[1], fv[2] });
    tri_fv.insert(tri_fv.end(), { fv[2], fv[3], fv[4] });
    tri_fv.insert(tri_fv.end(), { fv[4], fv[5], fv[0] });
    tri_fv.insert(tri_fv.end(), { fv[0], fv[2], fv[4] });
}

// 
void PLCc::get_triangles_naive(std::vector<uint32_t>& tri_fv) const {

    // DOES NOT WORK WITH SIMPLIFIRD CHAMFERED PLC

    uint32_t n_exact_vrts, n_LNC_vrts, n_BPT_vrts;

    size_t n_init_faces = faces.size();
    for (size_t i = 0; i < n_init_faces; i++) {
        std::vector<uint32_t> fv;
        get_face_vertices(faces[i], fv);

        n_exact_vrts = n_LNC_vrts = n_BPT_vrts = 0;

        for (uint32_t v : fv) {
            if (vertices[v]->isExplicit3D()) n_exact_vrts++;
            else if (vertices[v]->isLNC()) n_LNC_vrts++;
            else if (vertices[v]->isBPT()) n_BPT_vrts++;
#ifdef PLCC_DEBUG
            else {
                std::cout << "ERROR - triangulate_chamfered_plc: unexpected vertex type\n";
                assert(false);
            }
#endif
        }

        int32_t face_class = 100 * n_exact_vrts + 10 * n_LNC_vrts + n_BPT_vrts;
        switch (face_class) {

        case 300: // it is an input triangle.
            tri_fv.insert(tri_fv.end(), fv.begin(), fv.end());
            break;

        case 222: // obtained by chamfering 1 vertex from an input triangle
            get_trivrts_1av(fv, tri_fv);
            break;

        case 144: // obtained by chamfering 2 vertices from an input triangle
            get_trivrts_2av(fv, tri_fv);
            break;

        case 124: // obtained by chamfering 2 vertices and 1 edge from an input triangle
            get_trivrts_2av1ae(fv, tri_fv);
            break;

        case  66: // obtained by chamfering 3 vertices from an input triangle
            get_trivrts_3av(fv, tri_fv);
            break;

        case  46: // obtained by chamfering 3 vertices and 1 edge from an input triangle
            get_trivrts_3av1ae(fv, tri_fv);
            break;

        case  26: // obtained by chamfering 3 vertices and 2 edge from an input triangle
            get_trivrts_3av2ae(fv, tri_fv);
            break;

        case   6: // obtained by chamfering 3 vertices and 3 edge from an input triangle
            get_trivrts_3av3ae(fv, tri_fv);
            break;

#ifdef PLCC_DEBUG
        default:
            std::cout << "ERROR - triangulate_chamfered_plc: unexpected classification\n";
            assert(false);
#endif

        }

    }

}
