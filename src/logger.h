#include <chrono>

FILE* log_fp;
std::chrono::steady_clock::time_point time_point;

inline void startLogging(const char* fn) {
    
    char log_file_name[] = "delOpt_log.csv";
    char first_line[] = "Input_File, num_vrts, num_tets, "
                        "Elapsed_time(ms), Shortest_Edge, " 
                        "max_Energy_int, max_Energy_ext, "
                        "min_Face_Ang_int(DEG), max_Face_Ang_int(DEG), "
                        "min_Face_Ang_ext(DEG), max_Face_Ang_ext(DEG), "
                        "min_Dihed_Ang_int(DEG), max_Dihed_Ang_int(DEG), "
                        "min_Dihed_Ang_ext(DEG), max_Dihed_Ang_ext(DEG)"
                        "\n";

    if (fn != NULL) {
        log_fp = fopen(log_file_name, "r");
        if (log_fp == NULL) {
            log_fp = fopen(log_file_name, "w");
            fprintf(log_fp, "%s", first_line);
        }
        else {
            fclose(log_fp);
            log_fp = fopen(log_file_name, "a");
        }
        if (log_fp == NULL) ip_error("Can't open the file for logging!\n");

        size_t i;
        for (i = strlen(fn); i > 0; i--) if (fn[i - 1] == '\\' || fn[i - 1] == '/') break;
        fprintf(log_fp, "%s", fn + i);
    }
    else {
        log_fp = stdout;
        fprintf(log_fp, "%s", first_line);
    }

    time_point = std::chrono::steady_clock::now();
}

inline void logTimeChunk() {
    std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
    uint64_t ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_point).count();
    time_point = now;

    fprintf(log_fp, ", %zu", ms);
}

inline void logBoolean(bool b) {
    fprintf(log_fp, ", %s", b ? "True" : "False");
}

inline void logInteger(uint32_t n) {
    fprintf(log_fp, ", %u", n);
}

inline void logDouble(double d) {
    fprintf(log_fp, ", %g", d);
}

inline void logEmpty() {
    fprintf(log_fp, ", ");
}

inline void finishLogging() {
    fprintf(log_fp, "\n");
    if (log_fp != stdout) fclose(log_fp);
}

#ifdef _MSC_VER
#include <windows.h>
#include <psapi.h>

// To ensure correct resolution of symbols, add Psapi.lib to TARGETLIBS
// and compile with -DPSAPI_VERSION=1

double getPeakMegabytesUsed()
{
    HANDLE hProcess = OpenProcess(PROCESS_QUERY_INFORMATION | PROCESS_VM_READ, FALSE, GetCurrentProcessId());
    if (NULL == hProcess) return 0;

    PROCESS_MEMORY_COUNTERS pmc;
    double mem = 0;
    if (GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc)))
    {
        mem = pmc.PeakWorkingSetSize / 1048576.0;
    }

    CloseHandle(hProcess);
    return mem;
}

// Mem info in Mb
inline void logMemInfo()
{
    fprintf(log_fp, ", %.2f", getPeakMegabytesUsed());
}
#else

#include <sys/time.h>
#include <sys/resource.h>

inline void logMemInfo() {
    struct rusage r_usage;
    getrusage(RUSAGE_SELF, &r_usage);
    fprintf(log_fp, ", %.2f", r_usage.ru_maxrss); // bytes
}
#endif
