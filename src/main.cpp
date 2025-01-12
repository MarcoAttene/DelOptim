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

	public:
	std::vector<genericPoint*> out_vertices;
	std::vector<uint32_t> out_conn_vertices; // necessary to close the surface before computing the final CDT
	std::vector<std::vector<uint32_t>> out_faces;
	bool input_is_manifold;
	bool input_has_interior;

	chamfering_interface(inputPLC& _plc, double _epsilon, 
						 bool safe = false, bool simplify_cham_plc = true, 
						 bool print_surf = false, bool verbose = false) {
		cut_plc = new PLCc(_plc, _epsilon, safe, simplify_cham_plc, verbose);

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
		if(print_surf) cut_plc->saveTriangles(out_tri, "chamfered_plc.off"); // save triangulated faces

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
	void createSteinerCDT(inputPLC& plc, bool comp_min_PLC_dist, bool produce_output, bool verbose =false) {
		
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
		
		if (produce_output) mesh->saveConstrTrisToOFF("constrainedFaces.off", constr_tri_asCorners);

		if (comp_min_PLC_dist) mesh->compute_min_inputPLC_dist(constr_tri_asCorners); // minimum distance between cdt elemnts constrained on the PLC
	}

	inline bool is_defined() const { return mesh!=NULL; }
	inline double min_inputPLC_dist(){ return mesh->get_min_inputPLC_dist(); }
};

// 'plc' is a valid input PLC to the process. Validity is assumed but not verified!
// TetMesh* createSteinerCDT(inputPLC& plc, bool min_PLC_dist, bool produce_output, uint32_t& nct, bool verbose =false) {

// 	// Build a delaunay tetrahedrization of the vertices
// 	TetMesh* tin = new TetMesh;
// 	tin->init_vertices(plc.coordinates.data(), plc.numVertices());
// 	tin->tetrahedrize();

// 	// Build a structured PLC linked to the Delaunay tetrahedrization
// 	PLCx Steiner_plc(*tin, plc.triangle_vertices.data(), plc.numTriangles());
// 	Steiner_plc.segmentRecovery_HSi(!verbose);
// 	Steiner_plc.faceRecovery(!verbose);
// 	std::vector<bool> constr_tri_asCorners;
// 	Steiner_plc.markInnerTets_andGetConstrFaces(constr_tri_asCorners);
// 	if (produce_output) tin->saveConstrTrisToOFF("constrainedFaces.off", constr_tri_asCorners);
// 	nct = tin->countConstrTris(constr_tri_asCorners);
// 	if (min_PLC_dist) tin->compute_min_inputPLC_dist(constr_tri_asCorners); // minimum distance between cdt elemnts constrained on the PLC

// 	return tin;
// }

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

// ------------------ //
// Collect statistics //
// ------------------ //

uint64_t take_time(chrono_clock::time_point& time_point){
	chrono_clock::time_point now = chrono_clock::now();
	uint64_t time_laps = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_point).count();
	time_point = now;
	return time_laps;
}

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

int main(int argc, char* argv[])
{
	initFPU();

	char filename[2048];
#ifndef DEBUG
	if (argc < 2) {
		std::cout << "Mesher - Create a well-shaped tetrahedral mesh out of a triangulated OFF file.\n";
		std::cout << "USAGE: ./delmesher [-v][-l][-o][-d 8] filename.off\n";
		std::cout << "OPTIONS:\n";
		std::cout << "\t[-m max_vrt] -> set to max_vrt the maximum number of vertices of the mesh during the Delaunay refinement stage.\n";
		std::cout << "\t[-a] -> extract from the Delaunay refined mesh a vaild input and gives it to a CDT algorithm that produces a conforming and quasi-optimal volume mesh (output see below) \n";
		std::cout << "\t[-s] -> produce as output the constrained surface of the quasi-optimal-CDT (see below)\n";
		std::cout << "\t[-i] -> (forces -a activation) if the input encloses a volume, uses only internal points and constrained faces as input for the quasi-optimal-CDT.\n";
		std::cout << "\t[-v] -> verbose mode\n";
		std::cout << "\t[-l] -> logging mode\n";
		std::cout << "\t[-d exp (pos. integer)] -> avoid Delaunay Refinement if the min pt.s dist after chamfering is < 10^-exp\n";
		std::cout << "\t[-b] -> computes the lower bound for minimum distance between mesh elements (point - segment - triangle), otherwise only point-point distances are computed.\n";
		std::cout << "\t[-o] -> produces outputs related to Delaunay refinement algorithm (see below)\n";
		std::cout << "\t[-c] -> produces an intermediate output (see below)\n";
		std::cout << "OUTPUT:\n";
		std::cout << "\t when [-a] is activated produces a volumetric (DRCDT_mesh.tet) mesh.\n";
		std::cout << "\t when [-s] is activated produces a surface (constrainedFaces.off) meshes.\n";
		std::cout << "\t when [-o] is activated produces a volumetric (DR_mesh.tet) and a surface (DR_plcfaces.off) meshes.\n";
		std::cout << "\t when [-c] is activated produces a surface (chamfered_plc.off) mesh of the chamfered plc (intermediate construction).\n";
		std::cout << "RETURNS:\n";
		std::cout << "\t0 when the whole execution terminates correctly (also when an iperror occours)\n";
		std::cout << "\t10 when option -d is activated and min dist. is violated\n";
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

	for (int i = 1; i < argc; i++)
		if (argv[i][0] == '-') {
			if (argv[i][1] == 'd') { min_dist_exp = atoi(argv[++i]); continue; }
			if (argv[i][1] == 'm') { max_vrts = atoi(argv[++i]); continue; }
			for (int j = 1; j < strlen(argv[i]); j++) options += argv[i][j];
		}
		else memcpy(filename, argv[i], strlen(argv[i]) + 1);

	bool log_mode = (options.find('l') != std::string::npos);
	bool verbose_mode = (options.find('v') != std::string::npos);
	bool output_cdt = (options.find('a') != std::string::npos);
	bool only_interior_output_cdt = (options.find('i') != std::string::npos);
	if(only_interior_output_cdt) output_cdt = true;
	bool produce_outcdt_surf = (options.find('s') != std::string::npos);
	bool produce_cahm_surf = (options.find('c') != std::string::npos);
	bool produce_DR_output = (options.find('o') != std::string::npos);
	bool lowBnd_onMeshDist = (options.find('b') != std::string::npos);

	// Internal parameters ---
	double epsilon; // chamfering distance; set just before chamfering below.
	bool simplify_chamferd_plc = true; // DEFAULT: true. Try to remove uncessary edges while keeping non-acute angles
	bool safe_chamfering = false; // DEFAULT: false. If true use "safe" chamfering, but creates shorter edges
	double optim_ratio = 2.0; // DEFAULT: 2.0. Delaunay Refinement "tetrhadron shape"
	bool remove_slivers = true; // DEFAULT: true
	bool display_histogram = false; // If true prints (on the shell) angles histograms.
	char out4CDT_name[] = ""; // If not empty an input for a cdt is saved on an off file with this string name.
	bool comp_DelRef_CDT = true; // runs a CDT algorithm using the Delauny Refinement enriched skin.
	bool log_inputCDT_stats = true; // if log_mode is ON register stats of the CDT of the input PLC
	// ------------------------

	chrono_clock::time_point time_point = chrono_clock::now(); // start timing
	chrono_clock::time_point time_zero = time_point;

	if (log_mode) startLogging(filename);

	// Load a valid PLC from file
	inputPLC plc;
	plc.initFromFile(filename, verbose_mode);

	if (!verbose_mode) std::cout<<"\ninput_file: "<<filename<<"\n";
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
	chamfering_interface cham(plc, epsilon, safe_chamfering, simplify_chamferd_plc, produce_cahm_surf, verbose_mode);
	bool input_enclose_vol = cham.input_has_interior;
	if (log_mode) log_chamferPLC_stats(cham, time_point);

	TetMesh tin; // Dealunay Tetrahedrization of the set of vertices of the chamfered PLC
	tin.init_vertices(cham.out_vertices);

	if (verbose_mode) std::cout << "Delaunizing vertices...\n";
	double BBox_len = tin.addBoundingBoxVertices(); // Adds to tin.vertices 8 new vertices defining a bounding box, by default dist = 1.0
	tin.tetrahedrize();

	if (log_mode) advance_ProcessLogging("Del_vertices");

	Tetrahedrization mesh; // Delaunay Refinement data structure

	// Copy the DT to the Delaunay Refinement data structure
	double closest_dist = mesh.initFromVerticesAndTets(tin.vertices, tin.tet_node, BBox_len);
	// now closest_dist take into account only vertex-vertex diatnce
	if (verbose_mode) std::cout << "Distance of closest points relative to bb diagonal: " << closest_dist << "\n";

	uint64_t time_DRinit = 0;
	if(log_mode) time_DRinit = take_time(time_point);

	cdt_interface input_cdt;
	if (lowBnd_onMeshDist){ 
		// OPTIONAL: Compute the minum distance between any two mesh elements (vertices, edges, triangles)
		//			 A cdt is needed to make this computation efficient.
		input_cdt.createSteinerCDT(plc, true, false); // needed to compute lower bound on generic mesh elements distances
		double min_PLC_dist = input_cdt.min_inputPLC_dist() / (3.0 * BBox_len); // the lower bound is 1/3 * min_PLC_dist normalized wrt mesh bounding box diagonal
		closest_dist = min(closest_dist, min_PLC_dist);
		if (verbose_mode) std::cout << "Distance of closest elems relative to bb diagonal: " << closest_dist << "\n";
	}

	if (log_mode){ 
		logDouble( closest_dist );
		logInteger( take_time(time_point) );
	}

	// OPTIONAL: if min_dist_exp has been given (!=UINT32_MAX) as user defined parameter,
	// avoid Delaunay Refinement and exit, if the min closest_dist is < 10^-min_dist_exp.
	if (min_dist_exp != UINT32_MAX){
		double min_dist = std::pow(10.0, (double)min_dist_exp * (-1.0));
		if (closest_dist < min_dist){ 
			if(log_mode) finishLogging();
			std::cout<<"\nPROGRAM ABORTED: Closest points are too close ( < "<< min_dist <<")\n\n\n";
			exit(10);
		}
	}
	
	if (log_mode) advance_ProcessLogging("DelRef_init");

	// Add PLC faces - at this stage faces are just collections of input triangles
	if (verbose_mode) std::cout << "Adding PLC faces...\n";
	mesh.addPLCFaces(cham.out_faces);

	if (log_mode){ 
		time_DRinit += take_time(time_point);
		logDouble(time_DRinit);
		advance_ProcessLogging("add_PLCfaces");
	}

	// Delaunay refinement algorithm
	// DR.1 - Split all missing segments while no more remain.
	// This does not modify the triangles in PLC faces, but updates their edges with
	// the subedges being computed.
	if (verbose_mode) std::cout << "Recovering segments...\n";
	
	mesh.recoverAllSegments();

	if (log_mode){ 
		logInteger( take_time(time_point) );
		advance_ProcessLogging("rec_seg");
	}

	// DR.2 - Create a local 2D Delaunay triangulation for each PLC face
	if (verbose_mode) std::cout << "Delaunizing faces...\n";
	mesh.delaunizePLCFaces();
	if (log_mode) advance_ProcessLogging("Del_faces");

	// DR.3 - Face recovery
	if (verbose_mode) std::cout << "Recovering faces...\n";
	mesh.recoverAllFaces();

	if (log_mode){ 
		logInteger( take_time(time_point) );
		advance_ProcessLogging("rec_faces");
	}

	// DR.4 - Tet optimization
	if (verbose_mode) {
		if(max_vrts==UINT32_MAX) std::cout << "Optimizing tets...\n";
		else std::cout << "Optimizing tets (max number of vertices = "<<max_vrts<<") ...\n";
	}
	mesh.optimizeTets(optim_ratio, remove_slivers, true, max_vrts);

	if (log_mode){ 
		logInteger( take_time(time_point) );
		advance_ProcessLogging("tet_optim");
	}

	// DR.5 - Classify internal/external tets after optimization
	// In case the input edges has even number of incident input triangles,
	// i.e. the input surface defines an "interior",
	// we create a CDT (or use an existing one) of the input surface 
	// to decide wherever an "otimized mesh" tetrahedron
	// is internal(fully contained) or external(not fully contained) wrt the input surface.
	if(input_enclose_vol){
		if (verbose_mode) std::cout<<"Input encloses a volume\n";
		if (!input_cdt.is_defined()) input_cdt.createSteinerCDT(plc, false, false); 
		markInternalTets(mesh, input_cdt.mesh);
		if (log_mode) advance_ProcessLogging("IntExt_class");
	}
	else markAllTetAsInternal(mesh); 

	if(log_mode) logInteger( input_cdt.time_cdt );
	
	if(verbose_mode){ 
		std::cout << "Elapsed time (ms): " << take_time(time_zero) << "\n";
		uint64_t bmem = getPeakRSS();
		std::cout << "Peak memory RSS (byte): " << bmem << "\n";
	}
	if (log_mode)  log_DelRef_mesh_stats(mesh, input_enclose_vol);
	else mesh.printReport(input_enclose_vol, "DelRef Mesh");
	
	if(display_histogram) make_histogram(mesh, "DR");

	if(produce_DR_output){
		mesh.saveTET("DR_mesh.tet");
		mesh.saveOFFInterface("DR_plcfaces.off");
	}

	// Delaunay Refinement + CDT

	if( comp_DelRef_CDT || strlen(out4CDT_name) != 0 ){

		if (verbose_mode) std::cout << "\nClosing chamfered PLC.\n";

		// Collects all Delaunay Refinement vertices and 
		// a Delaunay Refined triangulation of the chamfered plc, 
		// to be used as input for cdt
		std::vector<double> cdt_vrts;
		std::vector<uint32_t> cdt_tris;
		get_vrts_and_tris_for_cdt(mesh, cham.out_conn_vertices, cdt_vrts, cdt_tris, out4CDT_name, verbose_mode);

		if( comp_DelRef_CDT ){

			if (verbose_mode) std::cout << "\nComputing Delaunay Refined CDT.\n";

			inputPLC qo_plc; 
			qo_plc.initFromVectors(cdt_vrts.data(), cdt_vrts.size()/3, cdt_tris.data(), cdt_tris.size()/3, false);
			cdt_interface DR_cdt; 
			DR_cdt.createSteinerCDT(qo_plc, false, produce_outcdt_surf);
			Tetrahedrization stat_mesh;
			stat_mesh.initFromVerticesAndTets(DR_cdt.mesh->vertices, DR_cdt.mesh->tet_node);
			markAllTetAsInternal( stat_mesh );

			if(log_mode){
				log_mesh_stats(stat_mesh, take_time(time_point), DR_cdt.num_constr_tris);
				advance_ProcessLogging("final_CDT");
			}
			else stat_mesh.printReport(false, "CDT of partially Delaunay refined mesh");
			if(display_histogram) make_histogram(mesh, "CDT");
			if(display_histogram) make_histogram(mesh, "DR+CDT");
		}
		else{ if(log_mode) log_empty_mesh_stats(); }

	}

	if(log_mode && input_cdt.is_defined() && log_inputCDT_stats){
		Tetrahedrization stat_INmesh;
		stat_INmesh.initFromVerticesAndTets(input_cdt.mesh->vertices, input_cdt.mesh->tet_node);
		markAllTetAsInternal( stat_INmesh ); 
		if(strlen(out4CDT_name) != 0 && !comp_DelRef_CDT) log_empty_mesh_stats(); 
		log_mesh_stats(stat_INmesh, input_cdt.time_cdt, input_cdt.num_constr_tris);
		advance_ProcessLogging("reg_inCDT_stats");
	}
	else{ if(log_mode) log_empty_mesh_stats();  }
	
	if(log_mode) finishLogging();
	std::cout << "Execution correctly COMPLETED.\n\n\n";

#endif

	return 0;
}
