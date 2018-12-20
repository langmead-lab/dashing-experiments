#include "bonsai/hll/hll.h"
#include "bonsai/hll/bf.h"
#include "bonsai/hll/mh.h"
#include "bonsai/hll/aesctr/aesctr.h"
#include "bonsai/kspp/ks.h"
#include <getopt.h>
#include <omp.h>
#include <mutex>

using std::size_t;
using namespace sketch;

static const std::vector<size_t> DEFSZ {
    // Auto-generated with
    //print(", ".join(map(str, ((1 << i) for i in range(10, 28)))))
    // There are 18 sizes here.
    1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304
};

static const std::vector<double> DEF_OLAP_FRAC {
    //print("    " + ", ".join(("%.2f" % (0.05 * (i + 1)) for i in range(20))))
    // There are 10 elements here.
    0.05, 0.10, 0.25, 0.5, 0.75, 0.9
    //0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00
};

static const std::vector<size_t> DEF_SET_SZ {
    // Auto-generated with
    // print(", ".join(map(str, ((1 << i) for i in range(14, 33, 3)))))
    // The default number of pairs of sets is 21, as we repeat
    //
    16384, 131072, 1048576, 8388608, 67108864, 536870912
};

// With default parameters, there are 21 * 20 * 18 = 7560 experiments, with 7560 * 3 = 22680 lines

struct comb_t {
    size_t setsz1, setsz2, sketchsz;
    double olap_frac;
};

void emit_header(std::FILE *fp, size_t niter) {
    std::fprintf(fp, "#Type\tsize1\tsize2\ttrue ji\tsketch size\tbfji\tbfji_orig\thllji\tmhji\n");
    std::fflush(fp);
}

struct result_line_t {
    static constexpr size_t NCOLS = 4; // BF, HLL
    std::array<double, NCOLS> __attribute__((aligned(32))) values;
    result_line_t() {std::memset(this, 0, sizeof(*this));}
    result_line_t(const result_line_t &o) = delete;
    result_line_t(result_line_t &&o) = default;
    result_line_t &operator=(const result_line_t &o) = delete;
    result_line_t(std::initializer_list<double> vals) {
        if(vals.size() != NCOLS) throw std::runtime_error("Wrong size for initializer list");
        std::copy(vals.begin(), vals.end(), values.begin());
    }
    void printf(std::FILE *fp) const {
        std::fprintf(fp, "bfji: %lf. bfsb1: %lf. hlji: %lf\n",
                     values[0], values[1], values[2]);
    }
    result_line_t &operator+=(const result_line_t &o) {
        for(unsigned i(0); i < NCOLS; ++i) {
            values[i] += o.values[i];
        }
        return *this;
    }
    result_line_t &operator/=(double value) {
        value = 1./value;
        for(unsigned i(0); i < NCOLS; values[i++] *= value);
        return *this;
    }
    result_line_t &operator-=(double value) {
        for(auto &el: values) el -= value;
        return *this;
    }
    double &operator[](size_t index) {
        return values[index];
    }
    double operator[](size_t index) const {
        return values[index];
    }
    result_line_t &operator=(double value) {
        std::fill(values.begin(), values.end(), value);
        return *this;
    }
};
void print_result(const std::vector<result_line_t> &lines, ks::string &buffer, double true_ji, size_t sz1, size_t sz2, size_t sketchsz) {
    result_line_t tmp;
    for(const auto &el: lines) tmp += el;
    tmp /= lines.size();
    tmp -= true_ji;
    buffer.sprintf("BIAS\t%zu\t%zu\t%le\t%zu\t%le\t%le\t%le\t%le\n", sz1, sz2, true_ji, sketchsz,
                   tmp[0], tmp[1], tmp[2], tmp[3]);
    tmp = 0;
    for(const auto &el: lines)
        for(size_t i(0); i < result_line_t::NCOLS; ++i)
            tmp[i] += std::abs(el[i] - true_ji);
    tmp /= lines.size();
    buffer.sprintf("ABSERROR\t%zu\t%zu\t%le\t%zu\t%le\t%le\t%le\t%le\n", sz1, sz2, true_ji, sketchsz,
                   tmp[0], tmp[1], tmp[2], tmp[3]);
    tmp = 0;
    for(const auto &el: lines) {
        for(size_t i(0); i < result_line_t::NCOLS; ++i) {
            double tmpd = (el[i] - true_ji);
            tmpd *= tmpd;
            tmp[i] += tmpd;
        }
    }
    tmp /= lines.size();
    buffer.sprintf("MSE\t%zu\t%zu\t%le\t%zu\t%le\t%le\t%le\t%le\n", sz1, sz2, true_ji, sketchsz,
                   tmp[0], tmp[1], tmp[2], tmp[3]);
}

void perform_calculations(size_t setsz1, size_t setsz2, size_t sketchsz, size_t niter, double olap_frac, ks::string &buffer, bool redo_already_done=true) {
    buffer.clear();
#ifndef MAGIC_SEED
#define MAGIC_SEED 777uL
#endif
    const uint64_t seed = MAGIC_SEED;
    bf::bf_t bf1a(std::log2(sketchsz) + 3, 1, seed);
    hll::hll_t ha(std::log2(sketchsz));
    bf::bf_t bf1b(std::log2(sketchsz) + 3, 1, seed);
    hll::hll_t hb(std::log2(sketchsz));
    sketch::mh::RangeMinHash<uint64_t> rmha(1<<int(std::log2(sketchsz) - 3)), rmhb(1<<int(std::log2(sketchsz) - 3));
    const size_t shared_size = std::min(std::min(setsz1, setsz2) * 0.9, olap_frac * std::max(setsz2, setsz1));
    const double true_ji = static_cast<double>(shared_size) / (setsz1 + setsz2 - shared_size);
    std::vector<result_line_t> results;
    results.reserve(niter);
    for(size_t iternum(0); iternum < niter; ++iternum) {
        const auto aesseed = seed + setsz1 * setsz2 + iternum;
        aes::AesCtr<uint64_t, 8> gena(aesseed);
        for(size_t i = 0; i < shared_size; ++i) {
            uint64_t v = gena();
            //uint64_t v = gena() & (~UINT64_C(3));
            //assert((v & 0x3u) == 0);
            ha.addh(v);
            hb.addh(v);
            bf1a.addh(v);
            bf1b.addh(v);
            rmha.addh(v);
            rmhb.addh(v);
        }
        for(size_t i = 0; i < setsz1 - shared_size; ++i) {
            //uint64_t v = (gena() & (~UINT64_C(1))) | 2u;
            //assert((v & 0x3u) == 2);
            uint64_t v = gena();
            bf1a.addh(v);
            ha.addh(v);
            rmha.addh(v);
        }
        for(size_t i = 0; i < setsz2 - shared_size; ++i) {
            //uint64_t v = (gena() & (~UINT64_C(2))) | 1u;
            //assert((v & 0x3u) == 1);
            uint64_t v = gena();
            bf1b.addh(v);
            hb.addh(v);
            rmhb.addh(v);
        }
        ha.sum(); hb.sum();
        gena.seed(aesseed);
        double bfji1 = bf1a.jaccard_index(bf1b);
        double hfji1 = ha.jaccard_index(hb);
        double bfsb1 = bf1a.setbit_jaccard_index(bf1b);
        double mhji1 = rmha.jaccard_index(rmhb);
        results.emplace_back(std::initializer_list<double>{bfji1, bfsb1, hfji1, mhji1});
        ha.clear(); hb.clear(); bf1a.clear(); bf1b.clear();
    }
    print_result(results, buffer, true_ji, setsz1, setsz2, sketchsz);
}

void usage(char **argv=nullptr) {
    std::fprintf(stderr, "Usage: %s <flags>\n"
                         "-f\tadd to olap_fracs. (Defaults to 0.05, ..., 0.95 if empty). Can be used repeatedly.\n"
                         "-p\tNumber of threads [1]\n"
                         "-s\tAppend to sizes in bytes [equivalent to print(\", \".join(map(str, ((1 << i) for i in range(10, 28)))) if none specified]\n"
                         "-o\tWrite results to file at path. Default: /dev/stdout\n"
                         "-n\tSet sizes [equivalent to print(\", \".join(map(str, ((1 << i) for i in range(14, 36, 3))))) if none specified\n"
                         "-h\tEmit help menu\n", (argv ? *argv: const_cast<char *>("cmp_bf_hllf")));
    std::exit(1);
}

int main(int argc, char *argv[]) {
    std::vector<size_t> sizes_in_bytes;
    std::vector<size_t> set_sizes;
    std::vector<double> olap_fracs;
    std::FILE *ofp = stdout;
    int c, niter = 10, nthreads = 1, redo = true, simply_count_combs = false;
    while((c = getopt(argc, argv, "f:p:o:i:n:s:Nrh")) >= 0) {
        switch(c) {
            case 'N': simply_count_combs = true; break;
            case 'r': redo = false; break;
            case 'f': olap_fracs.push_back(std::atof(optarg)); break;
            case 'i': niter = std::atoi(optarg); break;
            case 's':
                sizes_in_bytes.push_back(std::strtoull(optarg, nullptr, 10));
                break;
            case 'p': nthreads = std::atoi(optarg); break;
            case 'n':
                set_sizes.push_back(std::strtoull(optarg, nullptr, 10));
                break;
            case 'o': if((ofp = std::fopen(optarg, "w")) == nullptr) throw 1; break;
            case 'h': usage(argv);
        }
    }
    if(sizes_in_bytes.empty()) sizes_in_bytes = DEFSZ;
    if(set_sizes.empty()) set_sizes = DEF_SET_SZ;
    if(olap_fracs.empty()) olap_fracs = DEF_OLAP_FRAC;
    std::vector<comb_t> combs;
    if(simply_count_combs) {
        size_t i = (((set_sizes.size() + 1) * (set_sizes.size())) >> 1) * sizes_in_bytes.size() * olap_fracs.size();
        std::fprintf(stderr, "Number of combinations: %zu\n", i);
        std::exit(EXIT_SUCCESS);
    }
    std::mutex m;
    for(size_t i(0); i < set_sizes.size(); ++i)
        for(size_t j(i); j < set_sizes.size(); ++j)
            for(const auto sketchsz: sizes_in_bytes)
                for(const auto olap_frac: olap_fracs)
                    combs.emplace_back(comb_t{set_sizes[i], set_sizes[j], sketchsz, olap_frac});
    emit_header(ofp, niter);
    std::vector<ks::string> buffers(nthreads);
    const int fn = fileno(ofp);
    omp_set_num_threads(nthreads);
    #pragma omp parallel for schedule(dynamic)
    for(size_t i = 0; i < combs.size(); ++i) {
        ks::string &to_use = buffers[omp_get_thread_num()];
        comb_t cb = combs[i];
        perform_calculations(cb.setsz1, cb.setsz2, cb.sketchsz, niter, cb.olap_frac, to_use, redo);
        {
            std::lock_guard<std::mutex> lock(m);
            to_use.write(fn);
        }
        to_use.clear();
    }
}
