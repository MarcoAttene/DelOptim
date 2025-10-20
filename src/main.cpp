#ifdef _MSC_VER // Workaround for known bug on MSVC
#define _HAS_STD_BYTE 0  // https://developercommunity.visualstudio.com/t/error-c2872-byte-ambiguous-symbol/93889
#endif
 
// #define USE_TETGEN // If uncommented, tetgen is used instead of our
					  // refinement algorithm. Works only with MSVC and CLang

// #define ONLY_LFS // If uncomment, computes the local feature size of the 
					// input triangulation (using a CDT algorithm) and exit.

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
	cout << "Creates a well-shaped tetrahedral mesh out of a triangulated "
			"closed surface.\n"
		 << "The input is required to be manifold, oriented and\nwith no "
			"degenerate or self-intersecting triangles.\n"
		 << "The mesh is produced by Delaunay-refining a box around the "
			"chamfered input.\n"
		 << "INPUT: an OFF file (filename.off) containing a closed "
			"triangulated surface.\n"
		 << "USAGE: ./delmesher [-v][-l][-e 8] filename.off\n"
		 << "OPTIONS:\n"
		 << "[-a]\tcreate an enriched CDT out of the Delaunay-refined mesh\n"
				"\tto produce an exactly conformal volume mesh.\n"
		 << "[-b]\tuse the Delaunay-refined mesh to produce an OFF file which\n"
				"\tcan be used as input for a CDT algorithm (see OUTPUT).\n"
		 << "[-c]\tcomputes the minimum LFS of the input triangulation.\n"
		 << "[-d]\tenables sliver removal during Delaunay refinement.\n"
		 << "[-e exp]\t(exp positive integer) do not perform Delaunay\n"
				"\trefinement if the LFS lower bound is < 10^-exp.\n"
		 << "[-m mv]\t(mv positive integer) interrupts Delaunay refinement\n"
				"\tas soon as it inserts max_vrt vertices.\n"
		 << "[-h]\tdisplay angles histograms.\n"
		 << "[-l]\tlogging mode (append a line to 'delOpt_log.csv').\n"
		 << "[-v]\tverbose mode.\n"
		 << "[-u]\tsaves the input after chamfering (see OUTPUT).\n"
		 << "[-w]\tsaves the outer boundary of the Delaunay-refined mesh "
				 "(see OUTPUT).\n"
		 << "[-x]\tsaves the Delaunay-refined mesh (see OUTPUT).\n"
		 << "[-y]\t(needs [-a]) saves the outer boundary of the enriched CDT "
				 "(see OUTPUT).\n"
		 << "[-z]\t(needs [-a]) saves the enriched CDT (see OUTPUT).\n"
		 << "OUTPUT:\n"
		 << "\tuse [-b] to produce a surface mesh. ('filename'_rebuilt.off)\n"
		 << "\tuse [-u] to produce a surface mesh. (chamfered_plc.off)\n"
		 << "\tuse [-w] to produce a surface mesh. (DR_interface.off)\n"
		 << "\tuse [-y] to produce a surface mesh. "
			"(DRCDT_constrainedFaces.off)\n"
		 << "\tuse [-x] to produce a volumetric mesh. (DR_mesh.tet)\n"
		 << "\tuse [-z] to produce a volumetric mesh. (DRCDT_mesh.tet)\n"
		 << "RETURNS:\n"
		 << "\t0 when the whole execution terminates correctly.\n"
		 << "\t10 when option -e is activated and min dist. is violated\n"
		 << "\n";
	exit(0);
}


#ifdef ONLY_LFS
void comp_bboxdiag_and_lfs(inputPLC& plc, bool verbose, bool log){
	// By default the input is enclosed in a bounding box whose diagonal is 
	// (1 + 2 * bbox_factor) * diagonal_of_the_minimal_bounding_box.
	double BBox_len = plc.get_BBox_diag();

	// Compute the minum distance between any two mesh elements (vertices, 
	// edges, triangles): a cdt is needed to make this computation efficient.
	cdt_interface input_cdt;
	input_cdt.createSteinerCDT(plc, true);
	double min_PLC_dist = input_cdt.inputPLC_LFS() / BBox_len;
	if (verbose) cout << "input PLC LFS relative to bb diagonal: " 
						<< min_PLC_dist << "\n";
	if (log){ 
		logDouble("lfs", min_PLC_dist);
		logTimeChunk("T_inPLC_lfs(ms)"); 
	}
	if (log) advance_ProcessLogging("lfs");
}
#endif

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
	bool comp_DelRef_CDT = (options.find('a') != string::npos);
	bool comp_inLFS = 	   (options.find('c') != string::npos);
	bool remove_slivers =  (options.find('d') != string::npos);
	bool histo = 		   (options.find('h') != string::npos);
	bool log_mode = 	   (options.find('l') != string::npos);
	bool verbose_mode =    (options.find('v') != string::npos);
	bool cham_out = 	   (options.find('u') != string::npos);
	bool DR_skin = 		   (options.find('w') != string::npos);
	bool DR_outmesh = 	   (options.find('x') != string::npos);
	bool DRCDT_skin = 	   (options.find('y') != string::npos);
	bool DRCDT_outmesh =   (options.find('z') != string::npos);

	// Internal parameters ----------------------------------------------------
	double bbox_factor = 1.0; // Bounding box diagonal elongation factor.
	double epsilon; // Chamfering distance; depends on bounding box diagonal.
	bool cham_simpl = true; // DEFAULT: true. Try to remove uncessary edges 
							// while keeping non-acute angles.
	bool cham_safe = true; // DEFAULT: false. If true use "safe" chamfering (
						   // algorithm with theoretical guarantees of removig
						   // all acute angles),  but creates shorter edges.
	double optim_ratio = 2.0; // DEFAULT: 2.0. Affects Delaunay Refinement 
							  // tetrhadra shape.
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

#ifdef ONLY_LFS
	
	comp_bboxdiag_and_lfs(plc, verbose_mode, log_mode); 

#elif USE_TETGEN
	
	Tetrahedrization mesh;
	mesh.initWithTetgen(plc.numVertices(), plc.coordinates.data(), 
				plc.numTriangles(), plc.triangle_vertices.data(), true, false);
	cout << "\n";

#else

	// Chamfering of the input PLC (epsilon controls the chamfering distance)
	// epsilon = plc.bbDiag() / 1000.0;	
	epsilon = DBL_MAX; // use the as great as possible epsilon
	chamfering_interface cham(cham_safe, cham_simpl, cham_out);
	cham.perform_chamfering(plc, epsilon, verbose_mode);
	if (log_mode) log_chamferPLC_stats(cham, plc);

	// TMP start -- for chamfering test --------
	delRef_interface DR_tmp(plc, cham, verbose_mode, log_mode);
	DR_tmp.init(comp_inLFS, min_dist_exp);
	if(log_mode){ 
		logInteger("Peak Mem (byte)", (uint64_t)getPeakRSS());
		finishLogging();
	}
	std::cout << "Execution correctly COMPLETED.\n\n\n";
	return 0;
	// TMP end -- for chamfering test --------

	// Delaunay Refinement of the chamfered PLC
	delRef_interface DR(plc, cham, verbose_mode, log_mode);
	DR.init(comp_inLFS, min_dist_exp);
	DR.recover_segments();
	DR.recover_faces();
	DR.optimize_tetrahedra(optim_ratio, remove_slivers, max_vrts);
	DR.intExt_classification();
	DR.get_mesh_statistics(log_mode, histo, DR_skin, DR_outmesh);

	// Delaunay Refinement + CDT
	if( comp_DelRef_CDT || strlen(out4CDT_name) != 0 )
		get_CDT_of_DRmesh(DR, cham, out4CDT_name, histo, DRCDT_skin, 
						DRCDT_outmesh, comp_DelRef_CDT, verbose_mode, log_mode);

	// Get statistics or agnles histogram of the input-CDT
	DR.get_inputCDT_statistics(histo, log_inputCDT_stats);
	
#endif

	uint64_t bmem = getPeakRSS();

	if(log_mode){ 
		logInteger("Peak Mem (byte)", bmem);
		finishLogging();
	}

	if(verbose_mode){ 
		std::cout << "Elapsed time (ms): " << take_time(main_time) << "\n";
		std::cout << "Peak memory RSS (byte): " << bmem << "\n";
	}

	std::cout << "Execution correctly COMPLETED.\n\n\n";
	return 0;
}
