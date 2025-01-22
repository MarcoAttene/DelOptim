#include <chrono>

FILE *log_fp, *log_prog;
std::chrono::steady_clock::time_point time_point;

inline void startLogging(const char* fn) {

    char log_prog_file_name[] = "delOpt_log_completed_steps.txt";
    
    char log_file_name[] = "delOpt_log.csv";

    if (fn != NULL) {

        // open a .csv file to collect statistics
        log_fp = fopen(log_file_name, "r");
        if (log_fp == NULL) { log_fp = fopen(log_file_name, "w"); }
        else {
            fclose(log_fp);
            log_fp = fopen(log_file_name, "a");
        }
        if (log_fp == NULL) ip_error("Can't open the file for logging!\n");

        // open a .txt file to monitor execution progresses (usefull for models that do not convege)
        log_prog = fopen(log_prog_file_name, "a");
        if (log_prog == NULL) { log_prog = fopen(log_prog_file_name, "w"); }
        else {
            fclose(log_prog);
            log_prog = fopen(log_prog_file_name, "a");
        }
        if (log_prog == NULL) ip_error("Can't open the file for logging progresses!\n");

        // get input model name
        size_t i;
        for (i = strlen(fn); i > 0; i--) if (fn[i - 1] == '\\' || fn[i - 1] == '/') break;
        fprintf(log_fp, "\nInput_File, %s", fn + i); fflush(log_fp);
        fprintf(log_prog, "\n%s", fn + i); fflush(log_prog);
    }
    else {
        log_prog = stdout;
        log_fp = stdout;
    }

    time_point = std::chrono::steady_clock::now(); // set "time zero"
}

inline void logTimeChunk(const char* name) {
    std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
    uint64_t ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_point).count();
    time_point = now;
    fprintf(log_fp, ", %s, %zu", name, ms); fflush(log_fp);
}
inline void skipTimeChunk() { time_point = std::chrono::steady_clock::now(); }

inline void logBoolean(const char* name, bool b) { fprintf(log_fp, ", %s, %s", name, b ? "True" : "False"); fflush(log_fp);  }
inline void logInteger(const char* name, uint32_t n) { fprintf(log_fp, ", %s, %u", name, n); fflush(log_fp); }
inline void logInteger(const char* name1, const char* name2, uint32_t n) { fprintf(log_fp, ", %s%s, %u", name1, name2, n); fflush(log_fp); }
inline void logInteger(const char* name, uint64_t n) { fprintf(log_fp, ", %s, %zu", name, n); fflush(log_fp); }
inline void logInteger(const char* name1, const char* name2, uint64_t n) { fprintf(log_fp, ", %s%s, %zu", name1, name2, n); fflush(log_fp); }
inline void logDouble(const char* name, double d) { fprintf(log_fp, ", %s, %g", name, d); fflush(log_fp); }
inline void logDouble(const char* name1, const char* name2, double d) { fprintf(log_fp, ", %s%s, %g", name1, name2, d); fflush(log_fp); }
inline void advance_ProcessLogging(const char* stage){ fprintf(log_prog, ", %s", stage); fflush(log_prog); }

inline void finishLogging() {
    if (log_fp != stdout) fclose(log_fp);
    if (log_prog != stdout) fclose(log_prog);
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
    fflush(log_fp);
}
#endif
