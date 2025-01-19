#ifdef _MSC_VER // Workaround for known bug on MSVC
#define _HAS_STD_BYTE 0  // https://developercommunity.visualstudio.com/t/error-c2872-byte-ambiguous-symbol/93889
#endif
 
// #define USE_TETGEN // actuallt does not work 

// #define DISP_PROGRESS // to display on the shell progresses during Delaunay Refinement alg.

#include <iostream>
#include <fstream>
#include "cdt.h"
#include "inputPLC.h"
#ifdef USE_TETGEN
#include "tetgen_interface.cpp"
#endif
#include <chrono>
#include "logger.h"
#include "getRSS.h"
#include "optim_cdt.hpp"

using namespace std;

#ifdef ISVISUALSTUDIO
#ifndef NDEBUG
#define DEBUG
#endif
#endif

#include "cham.h"

typedef std::chrono::steady_clock chrono_clock;

uint64_t take_time(chrono_clock::time_point& time_point){
	chrono_clock::time_point now = chrono_clock::now();
	uint64_t time_laps = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_point).count();
	time_point = now;
	return time_laps;
}

// -------------------------- //
// -- CHAMFERING INTERFACE -- //
//--------------------------- //

// It has 4 OPTIONS (boolean):
// - safe -> enables the "theoretically safe" chamfering strategy, which in turn creates shorter edges.
// - simplify_cham_plc -> enables merging unnecessary output edges as long as no acute angles appears. 
// - print_surf -> enables saving the triangulated output on the file chamfered_plc.off
// - verbose -> enables verbose mode
class chamfering_interface{

	private:
	PLCc* cut_plc;
	std::vector<uint32_t> out_tri; // to save triangles of the triangulated chamfered-PLC (as consecutive triple of indices wrt cut_plc->vertices)
	std::vector<uint32_t> uv; // Used Vertex 
	bool safe;
	bool simplify;
	bool save_output;

	public:
	std::vector<genericPoint*> out_vertices;
	std::vector<uint32_t> out_conn_vertices; // necessary to close the surface before computing the final CDT
	std::vector<std::vector<uint32_t>> out_faces;
	bool input_is_manifold;
	bool input_has_interior;

	chamfering_interface(bool _safe = false, bool _simplify = true, bool _save_output_surf = false) : safe(_safe), simplify(_simplify), save_output(_save_output_surf) {}

	void perform_chamfering(inputPLC& _plc, double _epsilon, bool verbose = false) {
		cut_plc = new PLCc(_plc, _epsilon, safe, simplify, verbose);

		cut_plc->get_triangles(out_tri); // computes a trinagulation of cut_plc (without modifing it) 
		if (verbose) std::cout << "Chamfered PLC triangulation COMPLETED\n\n";

		// uv is used to mark all vertices that have to be collected
		// becuse necessary to next steps of the algorithm.
		uv.resize(cut_plc->vertices.size(), UINT32_MAX);
		mark_input_vertices(); // needed to close the chamfered surface when Delaunay refinement is combined with CDT
		mark_out_tri_vertices(); // needed by Delaunay refinement

		uint32_t idx = 0;
		for (size_t i = 0; i < uv.size(); i++) if (uv[i] != UINT32_MAX) {
			uv[i] = idx++; // now uv stores new indexing
			out_vertices.push_back(cut_plc->vertices[i]);
			out_conn_vertices.push_back(cut_plc->ref_exp3D_vrt[i]); 
			// explicit points are all used, they are stored at the beginning of 
			// the vertices vector, so they do not change indexing.
		}

		out_faces.resize(out_tri.size() / 3);	
		for (size_t i = 0; i < out_tri.size() / 3; i++) {
			const uint32_t* otv = out_tri.data() + i*3;
			out_faces[i].assign({uv[ *otv ], uv[ *(otv+1) ], uv[ *(otv+2) ]});
		}

		// Produce some intermediate output mesh and save .off files 
		// cut_plc->save_rebuilded_input_after_chamfering(out_tri, "all_tris_chamf.off"); // Save chamfered input + complementar triangles to rebuild the input PLC
		// cut_plc->saveFaces(); // save polygonal faces
		if(save_output) cut_plc->saveTriangles(out_tri, "chamfered_plc.off"); // save triangulated faces

		input_is_manifold = cut_plc->input_plc_is_manifold();
		input_has_interior = cut_plc->input_plc_defines_interior();
	}

	inline void mark_input_vertices(){ for (size_t i = 0; i < cut_plc->n_in_vrts; i++) uv[i] = 1; }
	inline void mark_out_tri_vertices(){ for (uint32_t vi : out_tri) uv[ vi ] = 1; }

	inline uint32_t get_num_out_vrts(){ return (uint32_t)out_vertices.size(); }
	inline uint32_t get_num_out_faces(){ return (uint32_t)out_faces.size(); }
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
	
	cdt_interface() : mesh(NULL), time_cdt(0) {}

	// 'plc' is a valid input PLC to the process. Validity is assumed but not verified!
	void createSteinerCDT(inputPLC& plc, bool comp_min_PLC_dist, bool verbose =false) {
		
		chrono_clock::time_point time_zero = chrono_clock::now(); // start timing

		// Build a delaunay tetrahedrization of the vertices
		mesh = new TetMesh;
		mesh->init_vertices(plc.coordinates.data(), plc.numVertices());
		mesh->tetrahedrize();

		// Build a structured PLC linked to the Delaunay tetrahedrization
		PLCx Steiner_plc(*mesh, plc.triangle_vertices.data(), plc.numTriangles());
		Steiner_plc.segmentRecovery_HSi(!verbose);
		Steiner_plc.faceRecovery(!verbose);
		Steiner_plc.markInnerTets_andGetConstrFaces(constr_tri_asCorners);

		chrono_clock::time_point now = chrono_clock::now();
		time_cdt = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_zero).count();

		num_constr_tris = mesh->countConstrTris(constr_tri_asCorners);

		if (comp_min_PLC_dist) mesh->compute_min_inputPLC_dist(constr_tri_asCorners); // minimum distance between cdt elemnts constrained on the PLC
	}

	inline bool is_defined() const { return mesh!=NULL; }
	inline double min_inputPLC_dist(){ return mesh->get_min_inputPLC_dist(); }
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
	const pointType* v[4] = { t->v0()->getPoint(), t->v1()->getPoint(), t->v2()->getPoint(), t->v3()->getPoint() };
	getTetBarycenter(v, ccc);
	explicitPoint3D p(ccc[0], ccc[1], ccc[2]);
	cdt->vertices.push_back(&p);
	static uint64_t tet = 0;
	tet = cdt->searchTetrahedron(tet, cdt->numVertices() - 1);
	cdt->vertices.pop_back();
	return (cdt->mark_tetrahedra[tet>>2] == DT_IN);
}

void markInternalTets(Tetrahedrization& mesh, TetMesh* cdt) {
	
	//for (Tetrahedron* t : mesh.tets()) t->is_internal = isTetInternal(t, cdt);

	// Smarter method - exploit adjacencies to make location in CDT faster
	for (Tetrahedron* t : mesh.tets()) t->unmark<0>();
	Tetrahedra todo;
	todo.reserve(mesh.tets().size());
	Tetrahedron *s, *t = mesh.tets().front();
	todo.push_back(t); t->mark<0>();

	while (!todo.empty()) {
		t = todo.back();
		t->is_internal = isTetInternal(t, cdt);
		todo.pop_back();
		s = t->t0(); if (s != NULL && !s->isMarked<0>()) { todo.push_back(s); s->mark<0>(); }
		s = t->t1(); if (s != NULL && !s->isMarked<0>()) { todo.push_back(s); s->mark<0>(); }
		s = t->t2(); if (s != NULL && !s->isMarked<0>()) { todo.push_back(s); s->mark<0>(); }
		s = t->t3(); if (s != NULL && !s->isMarked<0>()) { todo.push_back(s); s->mark<0>(); }
	}

	for (Tetrahedron* t : mesh.tets()) t->unmark<0>();
}

inline void markAllTetAsInternal(Tetrahedrization& m){ for(Tetrahedron* t : m.tets()) t->is_internal = true; }

// ----------------------------- //
// Delaunay refinement interface //
// ----------------------------- //

class delRef_interface {
	private:
	inputPLC& plc;
	TetMesh tin;
	chamfering_interface& cham;

	chrono_clock::time_point& time_point;
	bool verbose, log;

	public:
	Tetrahedrization mesh;
	cdt_interface input_cdt;
	double BBox_len;
	double closest_dist;

	delRef_interface(inputPLC& _plc, chamfering_interface& _cham, chrono_clock::time_point& _time_point, bool _verbose, bool _log) : 
					 plc(_plc), verbose(_verbose), log(_log), cham(_cham), time_point(_time_point){};

	void init(bool fullLowBnd, uint32_t min_dist_exp){
		tin.init_vertices(cham.out_vertices);

		if (verbose) std::cout << "Delaunizing vertices...\n";
		BBox_len = tin.addBoundingBoxVertices(); // Adds to tin.vertices 8 new vertices defining a bounding box, by default dist = 1.0
		tin.tetrahedrize();

		if (log) advance_ProcessLogging("Del_vertices");

		// Copy the DT to the Delaunay Refinement data structure
		closest_dist = mesh.initFromVerticesAndTets(tin.vertices, tin.tet_node, BBox_len);
		// now closest_dist take into account only vertex-vertex diatnce
		if (verbose) std::cout << "Distance of closest points relative to bb diagonal: " << closest_dist << "\n";

		uint64_t time_DRinit = 0;
		if(log) time_DRinit = take_time(time_point);

		if (fullLowBnd){ 
			// OPTIONAL: Compute the minum distance between any two mesh elements (vertices, edges, triangles)
			//			 A cdt is needed to make this computation efficient.
			input_cdt.createSteinerCDT(plc, true); // needed to compute lower bound on generic mesh elements distances
			double min_PLC_dist = input_cdt.min_inputPLC_dist() / (3.0 * BBox_len); // the lower bound is 1/3 * min_PLC_dist normalized wrt mesh bounding box diagonal
			closest_dist = min(closest_dist, min_PLC_dist);
			if (verbose) std::cout << "Distance of closest elems relative to bb diagonal: " << closest_dist << "\n";
		}

		if (log){ 
			logDouble( closest_dist );
			logInteger( take_time(time_point) );
		}

		// OPTIONAL: if min_dist_exp has been given (!=UINT32_MAX) as user defined parameter,
		// avoid Delaunay Refinement and exit, if the min closest_dist is < 10^-min_dist_exp.
		if (min_dist_exp != UINT32_MAX){
			double min_dist = std::pow(10.0, (double)min_dist_exp * (-1.0));
			if (closest_dist < min_dist){ 
				if(log) finishLogging();
				std::cout<<"\nPROGRAM ABORTED: Closest points are too close ( < "<< min_dist <<")\n\n\n";
				exit(10);
			}
		}
		
		if (log) advance_ProcessLogging("DelRef_init");

		// Add PLC faces - at this stage faces are just collections of input triangles
		if (verbose) std::cout << "Adding PLC faces...\n";
		mesh.addPLCFaces(cham.out_faces);

		if (log){ 
			time_DRinit += take_time(time_point);
			logDouble(time_DRinit);
			advance_ProcessLogging("add_PLCfaces");
		}
	}

	void recover_segments(){
		// pt.1 - Split all missing segments while no more remain.
		// This does not modify the triangles in PLC faces, but updates their edges with
		// the subedges being computed.
		if (verbose) std::cout << "Recovering segments...\n";
		mesh.recoverAllSegments();
		if (log){ 
			logInteger( take_time(time_point) );
			advance_ProcessLogging("rec_seg");
		}
		// pt.2 - Create a local 2D Delaunay triangulation for each PLC face
		if (verbose) std::cout << "Delaunizing faces...\n";
		mesh.delaunizePLCFaces();
		if (log) advance_ProcessLogging("Del_faces");
	}

	void recover_faces(){
		if (verbose) std::cout << "Recovering faces...\n";
		mesh.recoverAllFaces();
		if (log){ 
			logInteger( take_time(time_point) );
			advance_ProcessLogging("rec_faces");
		}
	}
	
	void optimize_tetrahedra(double optim_ratio, bool remove_slivers, uint32_t max_vrts = UINT32_MAX){
		if (verbose) {
			if(max_vrts==UINT32_MAX) std::cout << "Optimizing tets...\n";
			else std::cout << "Optimizing tets (max number of vertices = "<<max_vrts<<") ...\n";
		}
		mesh.optimizeTets(optim_ratio, remove_slivers, true, max_vrts);
		if (log){ 
			logInteger( take_time(time_point) );
			advance_ProcessLogging("tet_optim");
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

		if(log) logInteger( input_cdt.time_cdt );
	}

};

// ------------------------------ //
// Logging and Collect statistics //
// ------------------------------ //

inline void log_inputPLC_stats(inputPLC& plc, chrono_clock::time_point& time_point){
	logInteger(plc.numVertices());
	logInteger(plc.numTriangles());
	logInteger( take_time(time_point) );
	advance_ProcessLogging("load_input");
}

inline void log_chamferPLC_stats(chamfering_interface& cham, chrono_clock::time_point& time_point){
	logBoolean(cham.input_is_manifold);
	logBoolean(!cham.input_has_interior);
	logInteger(cham.get_num_out_vrts());
	logInteger(cham.get_num_out_faces());
	logInteger( take_time(time_point) );
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

inline void logAngleStats(Tetrahedrization& mesh, bool input_encloses_vol = false){
	double maxEneIN, maxEneEX;
	mesh.maxTetEnergy(maxEneIN, maxEneEX);
	logDouble(maxEneIN); 
	if(input_encloses_vol) logDouble(maxEneEX); else logEmpty();
	
	double minFAIN, maxFAIN, minFAEX, maxFAEX, minDAIN, maxDAIN, minDAEX, maxDAEX;
	mesh.minMaxTetAngle(minFAIN, maxFAIN, minFAEX, maxFAEX, minDAIN, maxDAIN, minDAEX, maxDAEX);
	logDouble(minFAIN); 
	logDouble(maxFAIN); 
	if(input_encloses_vol) logDouble(minFAEX); else logEmpty();
	if(input_encloses_vol) logDouble(maxFAEX); else logEmpty();
	logDouble(minDAIN); 
	logDouble(maxDAIN); 
	if(input_encloses_vol) logDouble(minDAEX); else logEmpty();
	if(input_encloses_vol) logDouble(maxDAEX); else logEmpty();

	double avMinPlan, avMaxPlan, avMinDihed, avMaxDihed;
	mesh.averageAngles(avMinPlan, avMaxPlan, avMinDihed, avMaxDihed);
	logDouble(avMinPlan);
	logDouble(avMaxPlan);
	logDouble(avMinDihed);
	logDouble(avMaxDihed);
}

inline void log_DelRef_mesh_stats(Tetrahedrization& mesh, bool input_enclose_vol = false){
	logInteger( (uint32_t)mesh.num_vertices() );
	logInteger( (uint32_t)mesh.num_tetrahedra() );
	logInteger( (uint32_t)mesh.count_DelTris(true) ); // only constrained triangles
	logDouble( mesh.minEdgeLength() );
	logAngleStats(mesh, input_enclose_vol);
	advance_ProcessLogging("Optim");
}

inline void log_mesh_stats(Tetrahedrization& mesh, 
						   uint64_t time_ms,
						   uint32_t num_constr_tris,
						   bool input_enclose_vol = false){
	logInteger( time_ms );
	logInteger( (uint64_t) mesh.num_vertices() );
	logInteger( (uint64_t) mesh.num_tetrahedra() );
	logInteger( num_constr_tris );
	logAngleStats( mesh, input_enclose_vol );
}

inline void log_empty_mesh_stats(){ for(size_t i=0; i< 4 + 14; i++) logEmpty(); }

// -------------- //
// main functions //
// -------------- //

void optimizeChamferedInputPLC(){}

int main(int argc, char* argv[])
{
	initFPU();

	char filename[2048];
#ifndef DEBUG
	if (argc < 2) {
		std::cout << "Mesher - Create a well-shaped tetrahedral mesh out of a triangulated closed surface.\n";
		std::cout << "INPUT: an OFF file (filename.off) containing a closed triangulated surface.\n";
		std::cout << "USAGE: ./delmesher [-v][-l][-e 8] filename.off\n";
		std::cout << "OPTIONS:\n";
		std::cout << "[-a]\tafter the Delaunay refined mesh is created a CDT algorithm is used\n"
				  << "\tto produce a quasi-optimal volume mesh.\n";
		std::cout << "[-b]\textracts from the Delaunay refined mesh a triangulaed surface\n"
				  << "\tconforming the input surface, which can be used as input for\n"
				  << "\ta CDT algorithm (see OUTPUT).\n";
		std::cout << "[-c]\tcomputes the lower bound for minimum distance between\n"
				  << "\tmesh elements (point - segment - triangle),\n"
				  << "\totherwise only point-point distances is used.\n";
		std::cout << "[-d]\tenables sliver removal during Delaunay refinement stage.\n";
		std::cout << "[-e exp]\t(exp is positive integer) exit the program before Delaunay\n"
				  << "\trefinement if the minimum distance between mesh elements is < 10^-exp.\n";
		std::cout << "[-m mv]\t(mv is a positive integer) interrupts Delaunay refinament\n" 
				  << "\tas soon as the number of mesh vertices exceed max_vrt.\n";
		std::cout << "[-h]\t(terminal) display angles histograms.\n";
		std::cout << "[-l]\tlogging mode.\n";
		std::cout << "[-v]\tverbose mode.\n";
		std::cout << "[-u]\tenables chamfering output (see OUTPUT).\n";
		std::cout << "[-w]\tenables surface Delaunay refinement output (see OUTPUT).\n";
		std::cout << "[-x]\tenables volumetric Delaunay refinement output (see OUTPUT).\n";
		std::cout << "[-y]\t(needs [-a]) enables surface CDT of Delauny refinement output (see OUTPUT).\n";
		std::cout << "[-z]\t(needs [-a]) enables volumetric CDT of Delauny refinement output (see OUTPUT).\n";
		std::cout << "OUTPUT:\n";
		std::cout << "\t when [-b] is activated produces a surface mesh. ('filename'_rebuilt.off)\n";
		std::cout << "\t when [-u] is activated produces a surface mesh. (chamfered_plc.off)\n";
		std::cout << "\t when [-w] is activated produces a surface mesh. (DR_interface.off)\n";
		std::cout << "\t when [-y] is activated produces a surface mesh. (DRCDT_constrainedFaces.off)\n";
		std::cout << "\t when [-x] is activated produces a volumetric mesh. (DR_mesh.tet)\n";
		std::cout << "\t when [-z] is activated produces a volumetric mesh. (DRCDT_mesh.tet)\n";
		std::cout << "RETURNS:\n";
		std::cout << "\t0 when the whole execution terminates correctly (or when iperror occours)\n";
		std::cout << "\t10 when option -e is activated and min dist. is violated\n";
		std::cout << std::endl;
		return 0;
	}
	else strcpy(filename, argv[1]);
#else
	strcpy(filename, "..\\Input_file\\acute\\cup_fixed_fixed.off");
	//strcpy(filename, "..\\Input_file\\subdcube.off");
	//char filename[2048] = "..\\Input_file\\twocubes.off";

#endif

	// User interface options
	std::string options = "";
	uint32_t min_dist_exp = UINT32_MAX;
	uint32_t max_vrts = UINT32_MAX;
	char out4CDT_name[] = "";

	for (int i = 1; i < argc; i++)
		if (argv[i][0] == '-') {
			if (argv[i][1] == 'e') { min_dist_exp = atoi(argv[++i]); continue; }
			if (argv[i][1] == 'm') { max_vrts = atoi(argv[++i]); continue; }
			for (int j = 1; j < strlen(argv[i]); j++) options += argv[i][j];
		}
		else memcpy(filename, argv[i], strlen(argv[i]) + 1);

	if(options.find('b') !=  std::string::npos) strcpy(out4CDT_name, filename);
	bool comp_DelRef_CDT = 	 (options.find('a') != std::string::npos);
	bool fullLowBnd = 		 (options.find('c') != std::string::npos);
	bool remove_slivers = 	 (options.find('d') != std::string::npos);
	bool display_histogram = (options.find('h') != std::string::npos);
	bool log_mode = 		 (options.find('l') != std::string::npos);
	bool verbose_mode = 	 (options.find('v') != std::string::npos);
	bool cham_out = 		 (options.find('u') != std::string::npos);
	bool DR_skin = 			 (options.find('w') != std::string::npos);
	bool DR_outmesh = 		 (options.find('x') != std::string::npos);
	bool DRCDT_skin = 		 (options.find('y') != std::string::npos);
	bool DRCDT_outmesh = 	 (options.find('z') != std::string::npos);
	
	// Internal parameters ---
	double epsilon; // chamfering distance; set later beacause depends on bounding box diagonal.
	bool cham_simpl = true; // DEFAULT: true. Try to remove uncessary edges while keeping non-acute angles.
	bool cham_safe = false; // DEFAULT: false. If true use "safe" chamfering, but creates shorter edges.
	double optim_ratio = 2.0; // DEFAULT: 2.0. Delaunay Refinement "tetrhadron shape".
	bool log_inputCDT_stats = true; // if log_mode is ON register stats of the CDT of the input PLC.
	// ------------------------

	chrono_clock::time_point time_point = chrono_clock::now(); // start timing
	chrono_clock::time_point time_zero = time_point;

	if (log_mode) startLogging(filename);

	// Load a valid PLC from file
	inputPLC plc;
	plc.initFromFile(filename, verbose_mode);

	if (log_mode) log_inputPLC_stats(plc, time_point);

#ifdef USE_TETGEN
	Tetrahedrization mesh;
	mesh.initWithTetgen(plc.numVertices(), plc.coordinates.data(), plc.numTriangles(), plc.triangle_vertices.data(), true, false);
	for (Tetrahedron* t : mesh.tets()) t->is_internal = true;
	mesh.printReport();
	std::cout << std::endl;
#else

	// Chamfering of the input PLC
	epsilon = plc.bbDiag() / 1000.0;	// Use bounding-box-diagonal/1000 as chamfering distance
	chamfering_interface cham(cham_safe, cham_simpl, cham_out);
	cham.perform_chamfering(plc, epsilon, verbose_mode);
	bool input_enclose_vol = cham.input_has_interior;
	if (log_mode) log_chamferPLC_stats(cham, time_point);

	// Delaunay Refinement of the chamfered PLC
	delRef_interface DR(plc, cham, time_point, verbose_mode, log_mode);
	DR.init(fullLowBnd, min_dist_exp);
	DR.recover_segments();
	DR.recover_faces();
	DR.optimize_tetrahedra(optim_ratio, remove_slivers, max_vrts);
	DR.intExt_classification();
	
	if(verbose_mode){ 
		std::cout << "Elapsed time (ms): " << take_time(time_zero) << "\n";
		uint64_t bmem = getPeakRSS();
		std::cout << "Peak memory RSS (byte): " << bmem << "\n";
	}
	if (log_mode)  log_DelRef_mesh_stats(DR.mesh, input_enclose_vol);
	else DR.mesh.printReport(input_enclose_vol, "DelRef Mesh");
	
	if(display_histogram) make_histogram(DR.mesh, "DR");
	if(DR_skin) DR.mesh.saveOFFInterface("DR_interface.off");
	if(DR_outmesh) DR.mesh.saveTET("DR_mesh.tet");

	// Delaunay Refinement + CDT

	if( comp_DelRef_CDT || strlen(out4CDT_name) != 0 ){

		if (verbose_mode) std::cout << "\nClosing chamfered PLC.\n";

		// Collects all Delaunay Refinement vertices and 
		// a Delaunay Refined triangulation of the chamfered plc, 
		// to be used as input for cdt
		std::vector<double> cdt_vrts;
		std::vector<uint32_t> cdt_tris;
		// get_vrts_and_tris_for_cdt(mesh, cham.out_conn_vertices, cdt_vrts, cdt_tris, out4CDT_name, verbose_mode);
		get_vrts_and_tris_for_cdt(DR.mesh, cham.out_conn_vertices, cdt_vrts, cdt_tris, out4CDT_name, verbose_mode);

		if( comp_DelRef_CDT ){

			if (verbose_mode) std::cout << "\nComputing Delaunay Refined CDT.\n";

			inputPLC qo_plc; 
			qo_plc.initFromVectors(cdt_vrts.data(), cdt_vrts.size()/3, cdt_tris.data(), cdt_tris.size()/3, false);
			cdt_interface DR_cdt; 
			DR_cdt.createSteinerCDT(qo_plc, false);
			Tetrahedrization stat_mesh;
			stat_mesh.initFromVerticesAndTets(DR_cdt.mesh->vertices, DR_cdt.mesh->tet_node);
			markAllTetAsInternal( stat_mesh );

			if(log_mode){
				log_mesh_stats(stat_mesh, take_time(time_point), DR_cdt.num_constr_tris);
				advance_ProcessLogging("final_CDT");
			}
			else stat_mesh.printReport(false, "CDT of partially Delaunay refined mesh");
			if(display_histogram) make_histogram(stat_mesh, "DR+CDT");
			if(DRCDT_skin) DR_cdt.save_skin_toOFF("DRCDT");
			if(DRCDT_outmesh) DR_cdt.save_mesh_toTET("DRCDT");
		}
		else{ if(log_mode) log_empty_mesh_stats(); }

	}

	// if(log_mode && input_cdt.is_defined() && log_inputCDT_stats){
	if(log_mode && DR.input_cdt.is_defined() && log_inputCDT_stats){
		Tetrahedrization stat_INmesh;
		// stat_INmesh.initFromVerticesAndTets(input_cdt.mesh->vertices, input_cdt.mesh->tet_node);
		stat_INmesh.initFromVerticesAndTets(DR.input_cdt.mesh->vertices, DR.input_cdt.mesh->tet_node);
		markAllTetAsInternal( stat_INmesh ); 
		if(strlen(out4CDT_name) != 0 && !comp_DelRef_CDT) log_empty_mesh_stats(); 
		// log_mesh_stats(stat_INmesh, input_cdt.time_cdt, input_cdt.num_constr_tris);
		stat_INmesh.initFromVerticesAndTets(DR.input_cdt.mesh->vertices, DR.input_cdt.mesh->tet_node);
		advance_ProcessLogging("reg_inCDT_stats");
		if(display_histogram) make_histogram(stat_INmesh, "input CDT");
	}
	else{ if(log_mode) log_empty_mesh_stats();  }
	
	if(log_mode) finishLogging();
	std::cout << "Execution correctly COMPLETED.\n\n\n";

#endif

	return 0;
}
