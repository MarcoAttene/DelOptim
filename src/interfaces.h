#include <chrono>
#include "cdt.h"
#include "inputPLC.h"
#include "logger.h"
#include "getRSS.h"
#include "optim_cdt.hpp"
#include "cham.h"

typedef std::chrono::steady_clock chrono_clock;

uint64_t take_time(chrono_clock::time_point& time_point){
	chrono_clock::time_point now = chrono_clock::now();
	uint64_t time_laps = 
		std::chrono::duration_cast<
				std::chrono::milliseconds>(now - time_point).count();
	time_point = now;
	return time_laps;
}

// -------------------------- //
// -- CHAMFERING INTERFACE -- //
//--------------------------- //

// The first part of the algoritmh is called "chamfering" and consist on
// removing acute angles (both planar and solid) from the input triangulation 
// ('_plc').
// It produce a new plc ('cut_plc') where each facet is contained in one input 
// triangle. Output chamfered plc ('out_vertices' + 'out_faces') is then used 
// as input for the Delaunay refinement algorithm.
// 'out_conn_vertices' contains indices of input vertices that have become 
// isolated during chamfering but that are necessery to close rebuilt the
// input surface to compute the CDT of the refined mesh.

// It is influenced by 3 command line OPTIONS (boolean):
// - safe: enables the "theoretically safe" chamfering strategy, 
//		   which in turn creates shorter edges.
// - simplify: enables merging unnecessary output edges as long as 
//			   no acute angles appear. 
// - save_output: enables saving the triangulated output on the file 
//				  'chamfered_plc.off'.
class chamfering_interface{

private:
	PLCc* cut_plc;
	std::vector<uint32_t> out_tri; // Used to access triangles of the 
								   // (triangulated) chamfered-PLC (as triple 
								   // of indices wrt cut_plc->vertices)
	std::vector<uint32_t> uv; // Indices Used Vertices (wrt cut_plc->vertices),
							  // those that are incident to at least an output
							  // triangle.
	double chamPLC_part_lfs = DBL_MAX;

	// Command line options related parameters:
	bool safe;
	bool simplify;
	bool save_output;

public:
	std::vector<genericPoint*> out_vertices;
	std::vector<std::vector<uint32_t>> out_faces;
	bool input_is_manifold;
	bool input_has_interior;
	std::vector<uint32_t> out_conn_vertices; // necessary to close the surface 
											 // before computing the final CDT

	chamfering_interface( bool _safe = false, 
						  bool _simplify = true, 
						  bool _save_output = false ) : 
		safe(_safe), simplify(_simplify), save_output(_save_output), 
		cut_plc(nullptr), input_is_manifold(false), input_has_interior(false) {}

	void perform_chamfering(const inputPLC& _plc, 
							double _epsilon, bool verbose) {

		cut_plc = new PLCc(_plc, _epsilon, safe, simplify, verbose);
		cut_plc->get_triangles(out_tri); // Fills 'out_tri' with the 
										 //	trinagulated 'cut_plc' faces.
		if (verbose) std::cout << "Chamfered PLC triangulation COMPLETED\n\n";

		// 'uv' is used to mark all vertices that are involved in next steps of 
		// the algorithm. UINT32_MAX denotes an unused vertex.
		uv.resize(cut_plc->vertices.size(), UINT32_MAX);
		mark_input_vertices(); // Needed to close the chamfered surface when 
							   // Delaunay refinement is combined with CDT.
		mark_out_tri_vertices(); // Needed by Delaunay refinement.
		uint32_t idx = 0;
		for (size_t i = 0; i < uv.size(); i++) if (uv[i] != UINT32_MAX) {
			uv[i] = idx++; // now uv stores new indexing
			out_vertices.push_back(cut_plc->vertices[i]);
			out_conn_vertices.push_back(cut_plc->ref_exp3D_vrt[i]); 
			// Explicit points (which are all and only input '_plc' vertices) 
			// are all used; they are stored at the beginning of 
			// 'cut_plc->vertices' vector, so they do not change indexing.
			// Unused vertices may be created during chamfering simplification 
			// step.
		}

		// Each out face is a triangle denoted by a vector of 3 indices.
		out_faces.resize(out_tri.size() / 3);	
		for (size_t i = 0; i < out_tri.size() / 3; i++) {
			const uint32_t* otv = out_tri.data() + i*3;
			out_faces[i].assign({uv[ *otv ], uv[ *(otv+1) ], uv[ *(otv+2) ]});
		}

		// Save triangulated faces of chamfered PLC
		if(save_output) cut_plc->saveTriangles(out_tri, "chamfered_plc.off"); 
		// Save chamfered input + complementar triangles to rebuild input PLC
		// cut_plc->save_rebuilded_input_after_chamfering(
		//			out_tri, "all_tris_chamf.off"); 
		// cut_plc->saveFaces(); // Save polygonal faces
		
		input_is_manifold = cut_plc->input_plc_is_manifold();
		input_has_interior = cut_plc->input_plc_defines_interior();
	}

	inline void mark_input_vertices() { 
		for(size_t i = 0; i < cut_plc->n_in_vrts; ) uv[i++] = 1; 
	}
	inline void mark_out_tri_vertices() { for(uint32_t vi : out_tri) uv[vi]=1; }

	inline uint32_t get_num_out_vrts(){ return (uint32_t)out_vertices.size(); }
	inline uint32_t get_num_out_faces(){ return (uint32_t)out_faces.size(); }
	inline double get_chamPLC_part_lfs(){
		if(chamPLC_part_lfs==DBL_MAX) chamPLC_part_lfs = cut_plc->get_part_lfs(); 
		return chamPLC_part_lfs; 
	}
};

// ------------- //
// CDT interface //
// ------------- //

class cdt_interface{

	private:
	PLCx* Steiner_plc;
	std::vector<bool> constr_tri_asCorners;
	
	public:
	TetMesh* mesh;
	uint64_t time_cdt;
	uint32_t num_constr_tris;
	
	cdt_interface() : Steiner_plc(nullptr), mesh(nullptr), 
						time_cdt(0), num_constr_tris(0) { }

	// 'plc' is a valid input PLC. Validity is assumed but not verified!
	void createSteinerCDT(inputPLC& plc, bool comp_min_PLC_dist, 
												bool verbose =false) {
 		// start timing
		chrono_clock::time_point cdt_time_point = chrono_clock::now();

		// Build a delaunay tetrahedrization of the vertices
		mesh = new TetMesh;
		mesh->init_vertices(plc.coordinates.data(), plc.numVertices());
		mesh->tetrahedrize();
		// Build a PLC linked to the Delaunay tetrahedrization
		PLCx Steiner_plc(*mesh,plc.triangle_vertices.data(),plc.numTriangles());
		Steiner_plc.segmentRecovery_HSi(!verbose);
		Steiner_plc.faceRecovery(!verbose);
		Steiner_plc.markInnerTets_andGetConstrFaces(constr_tri_asCorners);
		time_cdt = take_time(cdt_time_point);
		num_constr_tris = mesh->countConstrTris(constr_tri_asCorners);
		if (comp_min_PLC_dist) { 
			mesh->set_inputPLC_LFS(plc.coordinates, plc.triangle_vertices,
								constr_tri_asCorners); 
		}
	}

	void createSteinerCDT(std::vector<genericPoint*>& v, 
							std::vector<uint32_t>& t, bool verbose =false) {
		// start timing
		chrono_clock::time_point cdt_time_point = chrono_clock::now();

		// Build a delaunay tetrahedrization of the vertices
		mesh = new TetMesh;

		// DEBUG : use explicit vertices
		// std::vector<const genericPoint*> ev; 
		// for(const genericPoint* vv : v){ 
		// 	double x,y,z; vv->getApproxXYZCoordinates(x,y,z);
		// 	ev.push_back(new explicitPoint3D(x,y,z));
		// }
		// mesh->vertices.insert(mesh->vertices.end(), ev.begin(), ev.end());

		mesh->vertices.insert(mesh->vertices.end(), v.begin(), v.end());
		mesh->inc_tet.resize(v.size(), UINT64_MAX);
		mesh->marked_vertex.resize(v.size(), 0);

		mesh->tetrahedrize();

		// Build a PLC linked to the Delaunay tetrahedrization
		PLCx Steiner_plc(*mesh, t.data(), t.size()/3);
		Steiner_plc.segmentRecovery_HSi(!verbose);
		Steiner_plc.faceRecovery(!verbose);
		Steiner_plc.markInnerTets_andGetConstrFaces(constr_tri_asCorners);

		time_cdt = take_time(cdt_time_point);

		num_constr_tris = mesh->countConstrTris(constr_tri_asCorners);
	}

	inline bool is_defined() const { return mesh!=NULL; }
	inline double inputPLC_LFS(){ return mesh->get_inputPLC_LFS(); }
	inline void save_skin_toOFF(const char* out_name){ 
		if(!is_defined()) return;
		char name[2048];  strcpy(name, out_name);
		strcat(name, "_constrainedFaces.off");
		mesh->saveConstrTrisToOFF(name, constr_tri_asCorners); 
	}
	inline void save_mesh_toTET(const char* out_name){
		if(!is_defined()) return;
		char name[2048];  strcpy(name, out_name);  
		strcat(name, "_mesh.tet");
		mesh->saveTET(name); 
	}
};

// ---------------------------------- //
// Internal / External classification //
// ---------------------------------- //

bool isTetInternal(Tetrahedron* t, TetMesh* cdt) {
	double ccc[3];
	const pointType* v[4] = { t->v0()->getPoint(), t->v1()->getPoint(), 
								t->v2()->getPoint(), t->v3()->getPoint() };
	getTetBarycenter(v, ccc);
	explicitPoint3D p(ccc[0], ccc[1], ccc[2]);
	cdt->vertices.push_back(&p);
	static uint64_t tet = 0;
	tet = cdt->searchTetrahedron(tet, cdt->numVertices() - 1);
	cdt->vertices.pop_back();
	return (cdt->mark_tetrahedra[tet>>2] == DT_IN);
}

void markInternalTets(Tetrahedrization& mesh, TetMesh* cdt) {
	
	// exploit adjacencies to make location in CDT faster
	for (Tetrahedron* t : mesh.tets()) t->unmark<0>();
	Tetrahedra todo;
	todo.reserve(mesh.tets().size());
	Tetrahedron* t = mesh.tets().front();
	todo.push_back(t); t->mark<0>();

	while (!todo.empty()) {
		t = todo.back();
		t->is_internal = isTetInternal(t, cdt);
		todo.pop_back();
		Tetrahedra incTet = {t->t0(), t->t1(), t->t2(), t->t3()};
		for(Tetrahedron* s : incTet) if (s != nullptr && !s->isMarked<0>()) { 
			todo.push_back(s); 
			s->mark<0>(); 
		} 
	}

	for (Tetrahedron* t : mesh.tets()) t->unmark<0>();
}

inline void markAllTetAsInternal(Tetrahedrization& m){ 
	for(Tetrahedron* t : m.tets()) t->is_internal = true; 
}

// ------------------------------ //
// Logging and Collect statistics //
// ------------------------------ //

inline void log_inputPLC_stats(inputPLC& plc){
	logInteger("Input_nV", plc.numVertices());
	logInteger("Input_nTri", plc.numTriangles());
	logTimeChunk("T_init(ms)");
	advance_ProcessLogging("load_input");
}

inline void log_chamferPLC_stats(chamfering_interface& cham){
	logBoolean("manifold", cham.input_is_manifold);
	logBoolean("open", !cham.input_has_interior);
	logInteger("cham_nVrts", cham.get_num_out_vrts());
	logInteger("cham_nTris", cham.get_num_out_faces());
	logDouble("cham_part_LFS", cham.get_chamPLC_part_lfs());
	logTimeChunk("T_cham(ms)");
	advance_ProcessLogging("chamfering");
}

void make_histogram(Tetrahedrization& mesh, const char* name){
	std::vector<uint32_t> fah, dah;
	mesh.makeFaceAngleHistogram(fah, 12);
	mesh.makeDihedralAngleHistogram(dah, 12);
	std::cout << name << " FACE ANGLE HISTOGRAM\n";
	for (auto x : fah) std::cout << x << "\t";
	std::cout << "\n";
	std::cout << name << " DIHEDRAL ANGLE HISTOGRAM\n";
	for (auto x : dah) std::cout << x << "\t";
	std::cout << "\n";
}

inline void logAngleStats(Tetrahedrization& mesh, const char* name, 
										bool input_encloses_vol = false){
	double maxEneIN, maxEneEX;
	mesh.maxTetEnergy(maxEneIN, maxEneEX);
	if(input_encloses_vol){ 
		logDouble(name,"_MaxEne_int",maxEneIN); 
		logDouble(name,"_MaxEne_ext",maxEneEX);
	}
	else logDouble(name,"_MaxEne",maxEneIN); 
	
	double minFAIN, maxFAIN, minFAEX, maxFAEX; 
	double minDAIN, maxDAIN, minDAEX, maxDAEX;
	mesh.minMaxTetAngle( minFAIN, maxFAIN, minFAEX, maxFAEX, 
							minDAIN, maxDAIN, minDAEX, maxDAEX );
	if(input_encloses_vol){ 
		logDouble(name,"_mFA_int(DEG)",minFAIN); 
		logDouble(name,"_MFA_int(DEG)",maxFAIN); 
		logDouble(name,"_mFA_ext(DEG)",minFAEX);
		logDouble(name,"_MFA_ext(DEG)",maxFAEX);
		logDouble(name,"_mDA_int(DEG)",minDAIN); 
		logDouble(name,"_MDA_int(DEG)",maxDAIN); 
		logDouble(name,"_mDA_ext(DEG)",minDAEX);
		logDouble(name,"_MDA_ext(DEG)",maxDAEX);
	}
	else{
		logDouble(name,"_mFA(DEG)",minFAIN); 
		logDouble(name,"_MFA(DEG)",maxFAIN); 
		logDouble(name,"_mDA(DEG)",minDAIN); 
		logDouble(name,"_MDA(DEG)",maxDAIN); 
	}

	double avMinPlan, avMaxPlan, avMinDihed, avMaxDihed;
	mesh.averageAngles(avMinPlan, avMaxPlan, avMinDihed, avMaxDihed);
	logDouble(name,"_av_mFA(DEG)",avMinPlan);
	logDouble(name,"_av_MFA(DEG)",avMaxPlan);
	logDouble(name,"_av_mDA(DEG)",avMinDihed);
	logDouble(name,"_av_MDA(DEG)",avMaxDihed);
}

inline void log_mesh_stats(Tetrahedrization& mesh, 
						   const char* mesh_name,
						   uint64_t time_ms,
						   uint32_t num_constr_tris,
						   bool input_enclose_vol = false){
	logInteger("T_", mesh_name, time_ms );
	logInteger(mesh_name, "_nVrts", (uint64_t) mesh.num_vertices() );
	logInteger(mesh_name, "_nTets", (uint64_t) mesh.num_tetrahedra() );
	logInteger(mesh_name, "_nConstrTris", num_constr_tris );
	logAngleStats( mesh, mesh_name, input_enclose_vol );
}

// ----------------------------- //
// Delaunay refinement interface //
// ----------------------------- //

class delRef_interface {
	private:
	inputPLC& plc;
	TetMesh tin;
	chamfering_interface& cham;

	bool verbose, log;

	public:
	Tetrahedrization mesh;
	cdt_interface input_cdt;
	double BBox_len;
	double closest_dist;

	delRef_interface(inputPLC& _plc, chamfering_interface& _cham, 
					bool _verbose, bool _log ) : 
					 plc(_plc), cham(_cham),
					 verbose(_verbose), log(_log),  
					 BBox_len(0), closest_dist(0) { };

	void init(double bbox_factor, bool comp_inLFS, uint32_t min_dist_exp){
		// start timing
		chrono_clock::time_point internal_time_point = chrono_clock::now(); 

		tin.init_vertices(cham.out_vertices);
		if (verbose) std::cout << "Delaunizing vertices...\n";

		// Adds to tin.vertices 8 new vertices defining a bounding box, 
		// by default dist = 1.0
        std::vector<double> bbox_coords;
        plc.getBoundingBoxVertices(bbox_coords, bbox_factor);
        BBox_len = plc.bbDiag() * (1 + 2*bbox_factor);
		tin.addBoundingBoxVertices(bbox_coords);
		tin.tetrahedrize();
		if (log) advance_ProcessLogging("Del_vertices");

		assert( mesh.checkConnectivity() == 0 ); // DEBUG
		assert( mesh.checkDelaunayCondition() ); // DEBUG

		// Copy the DT to the Delaunay Refinement data structure
		closest_dist = mesh.initFromVerticesAndTets(tin.vertices, tin.tet_node);
		// now closest_dist take into account only vertex-vertex diatnce
		if (verbose) std::cout << "Distance of tet-mesh closest points relative"
								  " to bb diagonal: " 
							   << closest_dist / BBox_len << "\n";
		if (verbose) std::cout << "Distance of chamfered closest points "
								  "relative to bb diagonal: " 
							<< cham.get_chamPLC_part_lfs() / BBox_len << "\n";

		uint64_t time_DRinit;
		if(log){ 
			// 1st partial initialization time
			time_DRinit = take_time(internal_time_point); 
			skipTimeChunk(); // sincronize logging clock
		}

        double min_PLC_dist = -1.0;
		if (comp_inLFS) { 
			// OPTIONAL: Compute the minum distance between any two mesh 
			//			 elements (vertices, edges, triangles)
			//			 A cdt is needed to make this computation efficient.
			input_cdt.createSteinerCDT(plc, true); 
			min_PLC_dist = input_cdt.inputPLC_LFS() / BBox_len;
			if (verbose) std::cout << "input PLC min-LFS relative to bb "
									  "diagonal: " << min_PLC_dist << "\n";
		}

		if (log){
            logDouble("INplc lfs", min_PLC_dist); 
			logDouble("closest_dist", closest_dist);
			logTimeChunk("T_closestDist(ms)"); // time to compute closest_dist
			take_time(internal_time_point); // sincronize internal clock
		}

		// OPTIONAL: if min_dist_exp has been given (!=UINT32_MAX) as user 
		//			 defined parameter, avoid Delaunay Refinement and exit, 
		//			 if the min closest_dist is < 10^-min_dist_exp.
		if (min_dist_exp != UINT32_MAX){
			double min_dist = std::pow(10.0, (double)min_dist_exp * (-1.0));
			if (closest_dist < min_dist){ 
				if(log) finishLogging();
				std::cout<<"\nPROGRAM ABORTED: lower bound on closest points "
						   "violated ( < "<< min_dist <<")\n\n\n";
				exit(10);
			}
		}
		
		if (log) advance_ProcessLogging("closest_dist");

		// Add PLC faces - faces are just collections of input triangles
		if (verbose) std::cout << "Adding PLC faces...\n";
		mesh.addPLCFaces(cham.out_faces);

		if (log){ 
			// 2nd partial initialization time
			time_DRinit += take_time(time_point); 
			logInteger("T_DRinit(ms)", time_DRinit);
			skipTimeChunk(); // sincronize logging clock
			advance_ProcessLogging("add_PLCfaces");
		}
	}

	void recover_segments(){
		// pt.1 - Split all missing segments while no more remain.
		// This does not modify the triangles in PLC faces, but updates their 
		// edges with the subedges being computed.
		if (verbose) std::cout << "Recovering segments...\n";
		mesh.recoverAllSegments();
		if (log){ 
			logTimeChunk( "T_DRrs(ms)" );
			advance_ProcessLogging("rec_seg");
		}

		assert( mesh.checkConnectivity() == 0 ); // DEBUG
		assert( mesh.checkDelaunayCondition() ); // DEBUG

		// pt.2 - Create a local 2D Delaunay triangulation for each PLC face
		if (verbose) std::cout << "Delaunizing faces...\n";
		mesh.delaunizePLCFaces();

		assert( mesh.checkConnectivity() == 0 ); // DEBUG
		assert( mesh.checkAllFaces() ); // DEBUG
		assert( mesh.checkDelaunayCondition() ); // DEBUG

		if (log){ 
			logTimeChunk( "T_DRdelf(ms)" );
			advance_ProcessLogging("Del_faces");
		}
	}

	void recover_faces(){
		if (verbose) std::cout << "Recovering faces...\n";
		mesh.recoverAllFaces();

		assert( mesh.checkConnectivity() == 0 ); // DEBUG
		assert( mesh.checkAllFaces() ); // DEBUG
		assert( mesh.compareDelTris(true) ); // DEBUG
		assert( mesh.checkDelaunayCondition() ); // DEBUG

		if (log){ 
			logTimeChunk("T_DRrf(ms)");
			advance_ProcessLogging("rec_faces");
		}
	}
	
	void optimize_tetrahedra(double optim_ratio, bool remove_slivers, uint32_t max_vrts = UINT32_MAX){
		if (verbose) {
			if(max_vrts==UINT32_MAX) std::cout << "Optimizing tets...\n";
			else std::cout << "Optimizing tets (max number of vertices = "<<max_vrts<<") ...\n";
		}
		mesh.optimizeTets(optim_ratio, remove_slivers, true, max_vrts);

		assert( mesh.checkConnectivity() == 0 ); // DEBUG
		assert( mesh.checkAllFaces() ); // DEBUG
		assert( mesh.compareDelTris(true) ); // DEBUG
		assert( mesh.checkDelaunayCondition() ); // DEBUG

		if (log){ 
			logTimeChunk("T_DRopt(ms)");
			advance_ProcessLogging("optim_tets");
		}
	}

	void intExt_classification(){
		// In case the input edges has even number of incident input triangles,
		// i.e. the input surface defines an "interior",
		// we create a CDT (or use an existing one) of the input surface 
		// to decide wherever an "otimized mesh" tetrahedron
		// is internal(fully contained) or external(not fully contained) wrt the input surface.
		if(cham.input_has_interior){
			if (verbose) std::cout<<"Input encloses a volume\n";
			if (!input_cdt.is_defined()) input_cdt.createSteinerCDT(plc, false); 
			markInternalTets(mesh, input_cdt.mesh);
			if (log) advance_ProcessLogging("IntExt_class");
		}
		else markAllTetAsInternal(mesh); 

		if(log) logTimeChunk("T_DRintExt(ms)");
	}

	void get_mesh_statistics(bool log, bool histo, bool DR_skin, bool DR_outmesh){
		bool input_enclose_vol = cham.input_has_interior;
		if (log){
			logInteger("DR_nVrts", (uint64_t)mesh.num_vertices() );
			logInteger("DR_nTets", (uint64_t)mesh.num_tetrahedra() );
			logInteger("DR_nConstrTris", (uint64_t)mesh.count_DelTris(true) ); // only constrained triangles
			logDouble("DR_Short_Edge", mesh.minEdgeLength() );
			logAngleStats(mesh, "DR", input_enclose_vol);
			advance_ProcessLogging("Optim");
		}
		else mesh.printReport(input_enclose_vol, "DelRef Mesh");
		if(histo) make_histogram(mesh, "DR");
		if(DR_skin) mesh.saveOFFInterface("DR_interface.off");
		if(DR_outmesh) mesh.saveTET("DR_mesh.tet");
	}

	void get_inputCDT_statistics(bool display_histogram, bool log){
		if(!display_histogram && !log) return;
		if(!input_cdt.is_defined()) return;
		Tetrahedrization stat_INmesh;
		stat_INmesh.initFromVerticesAndTets(input_cdt.mesh->vertices, input_cdt.mesh->tet_node);
		markAllTetAsInternal( stat_INmesh ); 
		if(log) {
			log_mesh_stats(stat_INmesh, "INCDT", input_cdt.time_cdt, input_cdt.num_constr_tris);
			advance_ProcessLogging("reg_inCDT_stats");
		}
		if(display_histogram) make_histogram(stat_INmesh, "input CDT");
	}

};

// ---------------------------- //
// CDT of Delaunay refined mesh //
// ---------------------------- //

void get_CDT_of_DRmesh( delRef_interface& DR, chamfering_interface& cham, 
						const char* out4CDT_name,
						bool display_histogram, 
						bool DRCDT_skin, bool DRCDT_outmesh,
						bool computeCDT, bool verbose, bool log){

	if (verbose) std::cout << "\nClosing chamfered PLC.\n";

	// start timing
	chrono_clock::time_point loc_time_point = chrono_clock::now();

	// Collects all Delaunay Refinement vertices and 
	// a Delaunay refined triangulation of the chamfered plc, 
	// to be used as input for cdt
	std::vector<genericPoint*> cdt_vrts;
	std::vector<uint32_t> cdt_tris;
	get_vrts_and_tris_for_cdt(DR.mesh, cham.out_conn_vertices, 
								cdt_vrts, cdt_tris, out4CDT_name, verbose);

	if( !computeCDT ) return;

	if (verbose) std::cout << "\nComputing Delaunay Refined CDT.\n";

	cdt_interface DR_cdt; 
	DR_cdt.createSteinerCDT(cdt_vrts, cdt_tris);
	Tetrahedrization stat_mesh;
	stat_mesh.initFromVerticesAndTets(DR_cdt.mesh->vertices, 
										DR_cdt.mesh->tet_node);
	markAllTetAsInternal( stat_mesh );

	if(log){
		uint64_t drcdt_time = take_time(loc_time_point);
		log_mesh_stats(stat_mesh, "DRCDT", drcdt_time, DR_cdt.num_constr_tris);
		advance_ProcessLogging("final_CDT");
	}
	else stat_mesh.printReport(false, "CDT of partially Delaunay refined mesh");
	
	if(display_histogram) make_histogram(stat_mesh, "DR+CDT");
	if(DRCDT_skin) DR_cdt.save_skin_toOFF("DRCDT");
	if(DRCDT_outmesh) DR_cdt.save_mesh_toTET("DRCDT");
}