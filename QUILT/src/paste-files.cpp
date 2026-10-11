#include <Rcpp.h>
#include <zlib.h>
#include <cstdio>
#include <string>
#include <vector>
#ifndef _WIN32
#include <sys/resource.h>
#endif

// read one line (without the trailing newline) from a gz (or plain text) file
// returns false if at end of file
static bool gz_read_line(gzFile f, std::string& line, std::vector<char>& buf) {
    line.clear();
    while (true) {
        if (gzgets(f, buf.data(), buf.size()) == Z_NULL) {
            return !line.empty();
        }
        line.append(buf.data());
        if (!line.empty() && line.back() == '\n') {
            line.pop_back();
            return true;
        }
    }
}

struct GzFiles {
    std::vector<gzFile> in;
    gzFile out_gz = NULL;
    FILE* out = NULL;
    ~GzFiles() {
        for (gzFile f : in) {
            if (f != NULL) {
                gzclose(f);
            }
        }
        if (out_gz != NULL) {
            gzclose(out_gz);
        }
        if (out != NULL) {
            fclose(out);
        }
    }
};

// [[Rcpp::export]]
void Rcpp_paste_line_files(
    const std::vector<std::string>& files,
    const std::string& output_file,
    const int nLines,
    const bool append,
    const bool gz_output
) {
    // paste files of nLines lines together line by line (tab separated)
    // input files can be gz compressed or plain text
    // output is either gz compressed (fast compression), or plain text
    GzFiles g;
    g.in.resize(files.size(), NULL);
    for (size_t i = 0; i < files.size(); i++) {
        g.in[i] = gzopen(files[i].c_str(), "rb");
        if (g.in[i] == NULL) {
            Rcpp::stop("Cannot open file " + files[i]);
        }
    }
    if (gz_output) {
        g.out_gz = gzopen(output_file.c_str(), append ? "ab1" : "wb1");
        if (g.out_gz == NULL) {
            Rcpp::stop("Cannot open file " + output_file);
        }
    } else {
        g.out = fopen(output_file.c_str(), append ? "ab" : "wb");
        if (g.out == NULL) {
            Rcpp::stop("Cannot open file " + output_file);
        }
    }
    std::vector<char> buf(65536);
    std::string line, out_line;
    for (int iLine = 0; iLine < nLines; iLine++) {
        out_line.clear();
        for (size_t i = 0; i < g.in.size(); i++) {
            if (!gz_read_line(g.in[i], line, buf)) {
                Rcpp::stop("Unexpected end of file " + files[i] + " when pasting files together");
            }
            if (i > 0) {
                out_line.push_back('\t');
            }
            out_line.append(line);
        }
        out_line.push_back('\n');
        if (gz_output) {
            if (gzwrite(g.out_gz, out_line.data(), out_line.size()) != (int) out_line.size()) {
                Rcpp::stop("Error writing to file " + output_file);
            }
        } else {
            if (fwrite(out_line.data(), 1, out_line.size(), g.out) != out_line.size()) {
                Rcpp::stop("Error writing to file " + output_file);
            }
        }
        if (iLine % 10000 == 0) {
            Rcpp::checkUserInterrupt();
        }
    }
    return;
}



// [[Rcpp::export]]
int Rcpp_get_max_open_files(const int wanted) {
    // try to raise the soft limit on the number of open files to wanted
    // (capped by the hard limit), and return the (new) soft limit
#ifdef _WIN32
    return 500;
#else
    struct rlimit rl;
    if (getrlimit(RLIMIT_NOFILE, &rl) != 0) {
        return 500;
    }
    if (rl.rlim_cur != RLIM_INFINITY && rl.rlim_cur < (rlim_t) wanted) {
        rlim_t new_cur = (rlim_t) wanted;
        if (rl.rlim_max != RLIM_INFINITY && new_cur > rl.rlim_max) {
            new_cur = rl.rlim_max;
        }
        struct rlimit rl_new = rl;
        rl_new.rlim_cur = new_cur;
        if (setrlimit(RLIMIT_NOFILE, &rl_new) == 0) {
            rl.rlim_cur = new_cur;
        }
    }
    if (rl.rlim_cur == RLIM_INFINITY || rl.rlim_cur > (rlim_t) wanted) {
        return wanted;
    }
    return (int) rl.rlim_cur;
#endif
}
