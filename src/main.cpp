#ifdef _MSC_VER // Workaround for known bug on MSVC
#define _HAS_STD_BYTE 0  // https://developercommunity.visualstudio.com/t/error-c2872-byte-ambiguous-symbol/93889
#endif
 
// #define USE_TETGEN // If uncommented, tetgen is used instead of our
					  // refinement algorithm. Works only with MSVC and Clang

// #define DISP_PROGRESS // If uncommented, progress during Delaunay Refinement
					  	 // segment recovery, face recovery and optimization
					  	 // phases are displaied on standard output. It is
					  	 // suggested to use -v (verbose) command line option.

#include <iostream>
#include <fstream>
#include "interfaces.h"

using namespace std;

#ifdef ISVISUALSTUDIO
#ifndef NDEBUG
#define DEBUG
#endif
#endif

#ifdef USE_TETGEN
#include "tetgen_interface.cpp"
#endif

void  printUsageAndExit() {
	cout << "Creates a well-shaped and conformal tetrahedral mesh out of a "
			"triangulated, input surface.\n"
		 << "The input is required to be manifold, oriented and with no "
			"degenerate\nor self-intersecting triangles.\n"
		 << "The mesh is produced by Delaunay-refining a bounding box around "
			"the chamfered input.\n"
		 << "INPUT: an OFF file (filename.off) containing a triangulated"
			" surface.\n"
		 << "USAGE: ./delmesher [-v][-l][-e 8] filename.off\n"
		 << "OPTIONS:\n"
		 << "[-a]\tavoid to perform last phase of the algorithm called "
		 	"'enriched CDT';\n\t the output Delaunay-refined mesh will have "
			"quality guarantees (min face angle > 14 deg),\n\t but the "
			"constrained face will not be conformal to input surface.\n"
		 << "[-b]\tsaves the input to the e'enriched CDT' phase.\n"
		 << "[-c]\tcomputes the minimum LFS of the input triangulation.\n"
		 << "[-d]\tenables sliver removal during Delaunay refinement.\n"
		 << "[-e ex]\t(ex positive integer) set a lower bound on LFS of the "
			"chmfered PLC,\n\texcluding from chamfering and refinement "
			"processes all input\n\tconstraints responsable for a LFS smaller"
			"than 10^-ex times the\n\tBBox diagonal. The more input constarints"
			" are ignored, the more\n\tthe output mesh will resemble a CDT of "
			"the input triangulation.\n"
		 << "[-m mv]\t(mv positive integer) interrupts Delaunay refinement as "
			"soon as\n\tmv vertices are inserted.\n"
		 << "[-h]\tdisplay angles histograms.\n"
		 << "[-l]\tlogging mode (append a line to 'delOpt_log.csv').\n"
		 << "[-v]\tverbose mode.\n"
		 << "[-u]\tsaves the input after chamfering (see OUTPUT).\n"
		 << "[-w]\tsaves the outer boundary of the intermesiate "
			"Delaunay-refined mesh (see OUTPUT).\n"
		 << "[-x]\tsaves the intermediate Delaunay-refined mesh (see OUTPUT).\n"
		 << "[-y]\t(needs [-a]) saves the outer boundary of output mesh "
			"(see OUTPUT).\n"
		 << "[-z]\t(needs [-a]) saves the output (see OUTPUT).\n"
		 << "OUTPUT:\n"
		 << "\tuse [-b] to produce a surface mesh. ('filename'_rebuilt.off)\n"
		 << "\tuse [-u] to produce a surface mesh. (chamfered_plc.off)\n"
		 << "\tuse [-w] to produce a surface mesh. (DR_interface.off)\n"
		 << "\tuse [-y] to produce a surface mesh. (constrainedFaces.off)\n"
		 << "\tuse [-x] to produce a volumetric mesh. (DR_mesh.tet)\n"
		 << "\tuse [-z] to produce a volumetric mesh. (out_mesh.tet)\n"
		 << "RETURNS:\n"
		 << "\t0 when the whole execution terminates correctly.\n"
		 << "\n";
	exit(0);
}

// void comp_bboxdiag_and_lfs(inputPLC& plc, bool verbose, bool log){
// 	// Compute the minum distance between any two mesh elements (vertices, 
// 	// edges, triangles): a cdt is needed to make this computation efficient.
// 	cdt_interface input_cdt;
// 	input_cdt.createSteinerCDT(plc, true, false);
// 	double l = input_cdt.inputPLC_LFS() / plc.get_BBox_diag();
// 	if (verbose) cout << "input PLC LFS relative to bb diagonal: " << l << "\n";
// 	if (log){ 
// 		logDouble("lfs", l);
// 		logTimeChunk("T_inPLC_lfs(ms)"); 
// 		advance_ProcessLogging("lfs");
// 	}
// }

void preparing_to_return(uint64_t time, uint64_t mem, bool log, bool verbose) {
	if(log){ 
		logInteger("Peak Mem (byte)", mem);
		finishLogging();
	}

	if(verbose) std::cout << "Elapsed time (ms): " << time << "\n"
							<< "Peak memory RSS (byte): " << mem << "\n";

	std::cout << "Execution correctly COMPLETED.\n\n\n";
}

// ---- //
// main //
// ---- //

int main(int argc, char* argv[])
{
	initFPU();

	// User interface options
	std::string options = "";
	uint32_t min_dist_exp = UINT32_MAX;
	uint32_t max_vrts = UINT32_MAX;
	char out4CDT_name[2048] = "";
	char filename[2048] = "";

#ifndef DEBUG
	if (argc < 2) printUsageAndExit();
#else
	strcpy(filename, "..\\Input_file\\acute\\cup_fixed_fixed.off");
	//strcpy(filename, "..\\Input_file\\subdcube.off");
	//char filename[2048] = "..\\Input_file\\twocubes.off";
#endif

	for (int i = 1; i < argc; i++)
		if (argv[i][0] == '-') {
			if (argv[i][1] == 'e') { min_dist_exp = atoi(argv[++i]); continue; }
			if (argv[i][1] == 'm') { max_vrts = atoi(argv[++i]); continue; }
			for (int j = 1; j < strlen(argv[i]); j++) options += argv[i][j];
		}
		else memcpy(filename, argv[i], strlen(argv[i]) + 1);

	if (strlen(filename) == 0) printUsageAndExit();

	if(options.find('b') != string::npos) strcpy(out4CDT_name, filename);
	bool enriched_CDT =   		(options.find('a') == string::npos);
	bool comp_inLFS = 	  		(options.find('c') != string::npos);
	bool remove_slivers = 		(options.find('d') != string::npos);
	bool histo = 		  		(options.find('h') != string::npos);
	bool log_mode = 	  		(options.find('l') != string::npos);
	bool verbose_mode =   		(options.find('v') != string::npos);
	bool cham_out = 	  		(options.find('u') != string::npos);
	bool DR_skin = 		  		(options.find('w') != string::npos);
	bool DR_outmesh = 	  		(options.find('x') != string::npos);
	bool enriched_CDT_skin = 	(options.find('y') != string::npos);
	bool enriched_CDT_outmesh = (options.find('z') != string::npos);

	// Internal parameters ----------------------------------------------------
	double bbox_factor = 1.0; // Bounding box diagonal elongation factor.
	double epsilon; // Chamfering distance; depends on bounding box diagonal.
	bool cham_simpl = true; // DEFAULT: true. Try to remove uncessary edges 
							// while keeping non-acute angles.
	bool cham_safe = false; // DEFAULT: false. If true uses "safe" chamfering
						   // algorithm (with theoretical guarantees of removig
						   // all acute angles, but creates shorter edges).
	double optim_ratio = 2.0; // DEFAULT: 2.0. Affects Delaunay Refinement 
							  // tetrhadra shape.
	double toll = DBL_MAX; // tollerance to identify small features.
	bool log_inputCDT_stats = (log_mode && true); // DEFAULT: true. If true 
								//register stats of the CDT of the input PLC.
	// ------------------------------------------------------------------------
	
	chrono_clock::time_point main_time = chrono_clock::now(); // start timing
	if (log_mode) startLogging(filename);

	// Load a valid PLC from file
	inputPLC plc;
	plc.initFromFile(filename, verbose_mode);
	plc.setBoundingBox(bbox_factor);

	if (log_mode) log_inputPLC_stats(plc);

#ifdef USE_TETGEN
	Tetrahedrization mesh;
	mesh.initWithTetgen( plc.numVertices(), (double*)plc.ptr_to_coordinates(), 
	 					 plc.numTriangles(), (uint32_t*)plc.ptr_to_tri_vrts(), true, true);
	for(Tetrahedron* t : mesh.tets()) t->is_internal = true;
	mesh.saveOFFInterface("tetgen_interface.off");
	
	if (log_mode){
		logInteger("tetgen_nVrts", (uint64_t)mesh.num_vertices() );
		logInteger("tetgen_nTets", (uint64_t)mesh.num_tetrahedra() );
		logDouble("tetgen_Short_Edge", mesh.minEdgeLength() );
		logAngleStats(mesh, "tetgen", false);
		advance_ProcessLogging("tetgen");
	}
	else mesh.printReport(false, "tetgen Mesh");
	if(histo) make_histogram(mesh, "tetgen");
#else

	cdt_interface input_cdt;
	bool empty_plc = false;
	if(comp_inLFS || min_dist_exp!=UINT32_MAX) {
		input_cdt.createSteinerCDT(plc, comp_inLFS, min_dist_exp!=UINT32_MAX);
	}

	if(min_dist_exp != UINT32_MAX) {
		toll = std::pow(10.0, -1.0 * min_dist_exp);
		input_cdt.mark_small_tris(toll, plc);
		empty_plc = (plc.count_ignored_tris() == plc.numTriangles());
			
		if(verbose_mode) std::cout<< plc.count_ignored_tris() << " input "
							<< "triangles ignored (LFS/BB_DIG < "<<toll<<").\n";
	}

	// Chamfering of the input PLC (epsilon controls the chamfering distance)
	// epsilon = plc.bbDiag() / 1000.0;	
	epsilon = DBL_MAX; // use the as great as possible epsilon
	chamfering_interface cham(cham_safe, cham_simpl, cham_out, empty_plc);
	cham.perform_chamfering(plc, epsilon, toll, verbose_mode);
	if (log_mode) log_chamferPLC_stats(cham, plc);

	if(min_dist_exp!=UINT32_MAX) {
		if(verbose_mode) std::cout<< plc.count_ignored_tris() << " input "
						<< "triangles ignored (LFS/BB_DIG < "<<toll<<").\n";

		// plc.save_triangles("used_constraints.off", 1);
		// plc.save_triangles("unused_constraints.off", 2);
	}
	if(log_mode) logInteger("Used Input tris", 
					(uint64_t)(plc.numTriangles() - plc.count_ignored_tris()));

	// Delaunay Refinement of the chamfered PLC
	delRef_interface DR(plc, cham, input_cdt, verbose_mode, log_mode);
	DR.init(comp_inLFS);
	DR.recover_segments();
	DR.recover_faces();
	DR.optimize_tetrahedra(optim_ratio, remove_slivers, max_vrts);
	DR.intExt_classification();
	DR.get_mesh_statistics(log_mode, histo, DR_skin, DR_outmesh);

	// Enriched CDT
	if( enriched_CDT || strlen(out4CDT_name) != 0 )
		get_enriched_CDT(DR, cham, plc, toll, out4CDT_name, histo, 
						 enriched_CDT_skin, enriched_CDT_outmesh, enriched_CDT, 
						 verbose_mode, log_mode);

	// Get statistics or angles histogram of the input-CDT
	DR.get_inputCDT_statistics(histo, log_inputCDT_stats);
	
#endif

	uint64_t bmem = (uint64_t)getPeakRSS();
	preparing_to_return( take_time(main_time), bmem, log_mode, verbose_mode);
	return 0;
}
