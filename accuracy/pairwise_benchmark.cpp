#include "bonsai/hll/hll.h"
#include "bonsai/hll/bf.h"
#include "bonsai/bonsai/include/setcmp.h"
#include "bonsai/bonsai/include/encoder.h"
#include <omp.h>
#include <unordered_map>
#include <set>
#include <new>
#include <unistd.h>

using namespace sketch;
using namespace hll;
using namespace bf;
using namespace bns;


#if 0
1. Run pairwise comparisons
    1. Python script for performing calculation
    2. Perhaps C++ executable
    3. Compare sketch sizes from 1024 to `1 << 24`
    4. Compare all 3 JI estimate methods from hll.
    5. Compare:
        1. Range of JIs
            1. Pairs already selected
        2. Range of ANIs
            1. Pairs need to be selected
#endif

struct skglob_t;
struct gpair_t {
    std::string a, b;
    bool nonempty() const {return a.size() && b.size();}
};

#ifndef MASHPATH
#define MASHPATH "/home-1/dbaker49@jhu.edu/work/db/code/dash/Mash/mash"
#endif
ks::string make_fname(const char *fname, unsigned hllp, unsigned k, const std::string &scratch_folder, const char *s="msh");

ks::string sketch_mash(const char *path, unsigned hllp, unsigned k, const std::string &scratch_folder);
ks::string sketch_bindash(const char *path, unsigned hllp, unsigned k, const std::string &scratch_folder);
double get_mash_ji(const char *path1, const char *path2);
double get_bindash_ji(const char *path1, const char *path2);

struct skglob_t {
    using hvec_t = std::vector<hll_t>;
    hvec_t         hll_;
    std::vector<std::string> mash_paths_;
    std::vector<std::string> bindash_paths_;
    khash_t(all) *hash_;
    const std::string &scratch_folder_;
    skglob_t(const skglob_t &) = delete;
    skglob_t& operator=(const skglob_t &) = delete;
    skglob_t(skglob_t &&o): hll_(std::move(o.hll_)), mash_paths_(std::move(o.mash_paths_)), bindash_paths_(std::move(o.bindash_paths_)), hash_(o.hash_), scratch_folder_(o.scratch_folder_)
    {
        o.hash_ = nullptr;
    }
    ~skglob_t() {
        kh_destroy(all, hash_);
    }
    skglob_t(const char *path, unsigned pstart, unsigned pend, unsigned k, const std::string &scratch_folder, bool bd_only):
        hash_(kh_init(all)), scratch_folder_(scratch_folder)
    {
        if(!bd_only) {
            std::generate_n(std::back_inserter(hll_), pend - pstart + 1, [&](){
                return hll_t(hll_.size() + pstart, ERTL_MLE, (JointEstimationMethod)ERTL_MLE);
            });
            LOG_DEBUG("Made %zu hlls\n", hll_.size());
            gzFile fp = gzopen(path, "rb");
            if(fp == nullptr) {
                std::fprintf(stderr, "Could not open file at %s\n", path);
                throw std::runtime_error("Could not open file at "s + path);
            }
            Encoder<> enc(k);
    
            int khr;
            LOG_DEBUG("About to encode\n");
            enc.for_each([&](u64 v) {
                kh_put(all, hash_, v, &khr);
                v = sketch::hll::hll_t::HashType()(v);
                for(auto &h: hll_) h.add(v);
            }, fp);
            gzclose(fp);
            LOG_DEBUG("Closed gzfile. Now making mash!\n");
            mash_paths_ = std::vector<std::string>(pend - pstart + 1);
        } // if !bd_only
        bindash_paths_ = std::vector<std::string>(pend - pstart + 1);
        LOG_DEBUG("mp size: %zu. bd paths: %zu\n", mash_paths_.size(), bindash_paths_.size());
        #pragma omp parallel for
        for(unsigned i = 0; i < bindash_paths_.size(); ++i) {
            if(!bd_only) {
                mash_paths_[i] = make_fname(path, i + pstart, k, scratch_folder).data();
                if(!isfile(mash_paths_[i])) {
                    sketch_mash(path, i + pstart, k, scratch_folder_).data();
                    LOG_DEBUG("Successfully made mash sketch at %s\n", mash_paths_[i].data());
                } else {
                    LOG_INFO("msh:Already created, not needed.\n");
                }
            }
            bindash_paths_[i] = make_fname(path, i + pstart, k, scratch_folder, "bdsh").data();
            sketch_bindash(path, i + pstart, k, scratch_folder_).data();
        }
        LOG_DEBUG("Finished generating sketches for path %s.\n", path);
    }
};


std::vector<gpair_t> get_pairs_from_file(const char *fn) {
    gzFile fp = gzopen(fn, "rb");
    char *line, *p, *q, *t, buf[1024];
    std::vector<gpair_t> ret;
    while((line = gzgets(fp, buf, sizeof(buf)))) {
        if(*line == '#') continue;
        std::fprintf(stderr, "Passing line is %s\n", line);
        p = line;
        while(std::isspace(*p)) ++p;
        t = p;
        q = p + 1;
        while(!std::isspace(*q)) ++q;
        *q = '\0';
        t = p;
        p = ++q;
        while(!std::isspace(*p)) ++p;
        *p = '\0';
        ret.emplace_back(gpair_t{t, q});
        assert(ret.back().nonempty());
    }
    gzclose(fp);
    return ret;
}


class FullAnalysis {
    unsigned k_, pstart_, pend_, ls_;
    std::vector<gpair_t> pairs_;
    std::set<std::string> pathset_;
    std::unordered_map<std::string, skglob_t> map_;
    const bool use_ji_for_buckets_;
    const std::string &scratch_folder_;
    const bool bd_only_;
public:
    unsigned nps() const {
        return pend_ - pstart_ + 1;
    }
    struct ResultBatch {
        std::vector<double> exact_;
        std::vector<std::vector<double>> hll_orig_dists_;
        std::vector<std::vector<double>> hll_mle_dists_;
        std::vector<std::vector<double>> hll_jmle_dists_;
        std::vector<std::vector<double>> mash_dists_;
        std::vector<std::vector<double>> bindash_dists_;
        ~ResultBatch();
        const std::vector<gpair_t> &paths_;
        unsigned pstart_, k_;
        const bool bdo_;
        ResultBatch(size_t npairs, unsigned number_sizes, const std::vector<gpair_t> &paths, unsigned pstart, unsigned k, bool bdonly): exact_(npairs), paths_(paths), pstart_(pstart), k_(k), bdo_(bdonly) {
            bindash_dists_.reserve(npairs);
            while(bindash_dists_.size() < npairs)
                bindash_dists_.emplace_back(number_sizes);
            if(!bdo_) {
                hll_orig_dists_.reserve(npairs);
                hll_mle_dists_.reserve(npairs);
                hll_jmle_dists_.reserve(npairs);
                mash_dists_.reserve(npairs);
                while(hll_orig_dists_.size() < npairs) {
                    hll_orig_dists_.emplace_back(number_sizes);
                    hll_mle_dists_.emplace_back(number_sizes);
                    hll_jmle_dists_.emplace_back(number_sizes);
                    mash_dists_.emplace_back(number_sizes);
                }
            }
            LOG_INFO("Allocated vectors: exact size %zu\n", exact_.size());
        }
        void write(std::FILE *fp) const {
            for(size_t i(0); i < paths_.size(); ++i)
                for(size_t j(0); j < bindash_dists_[0].size(); ++j)
                    bdo_ ? std::fprintf(fp, "%s\t%s\t%u\t%u\t%le\n", paths_[i].a.data(), paths_[i].b.data(), k_, unsigned(j) + pstart_, bindash_dists_[i][j])
                         : std::fprintf(fp, "%s\t%s\t%u\t%u\t%le\t%le\t%le\t%le\t%le\t%le\n", paths_[i].a.data(), paths_[i].b.data(), k_, unsigned(j) + pstart_, exact_[i], hll_orig_dists_[i][j], hll_mle_dists_[i][j], hll_jmle_dists_[i][j], mash_dists_[i][j], bindash_dists_[i][j]);
            std::fflush(fp);
        }
    };
    ~FullAnalysis();
    FullAnalysis(unsigned k, const char *pairpath, unsigned pstart, unsigned pend, unsigned linspace, bool use_ji=true, const std::string &scratch_folder="", bool bindashonly=false):
        k_(k), pstart_(pstart), pend_(pend), ls_(linspace), pairs_(get_pairs_from_file(pairpath)), use_ji_for_buckets_(use_ji), scratch_folder_(scratch_folder), bd_only_(bindashonly)
    {
        LOG_INFO("Made full_analysis with %zu pairs\n", pairs_.size());
    }
    void fill(unsigned k=0) {
    }
    ResultBatch run(unsigned k=0) {
        if(pathset_.empty()) {
            for(const auto &pair: pairs_)
                pathset_.insert(pair.a), pathset_.insert(pair.b);
        }
        if(pend_ < pstart_) throw std::runtime_error("pend must be >= pstart");
        if(pstart_ < 8) throw std::runtime_error("hllstart must be 8 or greater");
        if(pend_ > 32) throw std::runtime_error("hllend must be 32 or less");
        if(k) k_ = k;
        
        for(const auto &s: pathset_) {std::fprintf(stderr, "path: %s\n", s.data());}
        for(const auto &s: pathset_) {
            map_.emplace(s, skglob_t(s.data(), pstart_, pend_, k, scratch_folder_, bd_only_));
            std::fprintf(stderr, "map is size %zu of %zu\n", map_.size(), pathset_.size());
        }
        LOG_INFO("For k = %u, filled sketches with a total number of kv pairs in map_ of %zu\n", k, map_.size());
        ResultBatch ret(pairs_.size(), nps(), pairs_, pstart_, k, bd_only_);
        if(bd_only_) {
            #pragma omp parallel for
            for(unsigned i = 0; i < pairs_.size(); ++i) {
                auto &pair = pairs_[i];
                auto &f1 = map_.at(pair.a);
                auto &f2 = map_.at(pair.b);
                for(unsigned j = 0; j < ret.bindash_dists_[0].size(); ++j) {
                    ret.bindash_dists_[i][j] = get_bindash_ji(f1.bindash_paths_[j].data(), f2.bindash_paths_[j].data());
                }
            }
            goto cleanup;
        }
        #pragma omp parallel for
        for(unsigned i = 0; i < pairs_.size(); ++i) {
            LOG_INFO("Processing pair at index %u with tid = %i and paths %s, %s\n", i, omp_get_thread_num(), pairs_[i].a.data(), pairs_[i].b.data());
            auto &pair = pairs_[i];
            auto &f1 = map_.at(pair.a);
            auto &f2 = map_.at(pair.b);
            assert(f1.hash_);
            assert(f2.hash_);
            ret.exact_[i] = bns::jaccard_index(f1.hash_, f2.hash_);
            for(unsigned j = 0; j < ret.hll_orig_dists_[i].size(); ++j) {
                auto &h1 = f1.hll_[j], &h2 = f2.hll_[j];
                h1.set_estim(ORIGINAL);
                h1.set_jestim(ORIGINAL);
                h1.sum(); h2.sum();
                ret.hll_orig_dists_[i][j] = h1.jaccard_index(h2);
                h1.set_estim(ERTL_MLE);
                h1.set_jestim(ERTL_MLE);
                h1.sum(); h2.sum();
                ret.hll_mle_dists_[i][j] = h1.jaccard_index(h2);
                h1.set_estim(ERTL_JOINT_MLE);
                h1.set_jestim(ERTL_JOINT_MLE);
                h1.sum(); h2.sum();
                ret.hll_jmle_dists_[i][j] = h1.jaccard_index(h2);
                ret.mash_dists_[i][j] = get_mash_ji(f1.mash_paths_[j].data(), f2.mash_paths_[j].data());
                if(f1.bindash_paths_.at(j).empty() || f2.bindash_paths_.at(j).empty()) {
                    // std::fprintf(stderr, "Sizes of these paths are: %zu/%zu\n", f1.bindash_paths_[j].size(), f2.bindash_paths_[j].size());
                }
                ret.bindash_dists_[i][j] = get_bindash_ji(f1.bindash_paths_[j].data(), f2.bindash_paths_[j].data());
            }
        }
        cleanup:
        LOG_INFO("Calculated all distances for k = %u\n", k);
        map_.clear();
        return std::move(ret);
    }
    template<typename T, typename Functor>
    void krange_for_each(T it1, T it2, const Functor &func) {
        static_assert(std::is_integral<std::decay_t<decltype(*it1)>>::value, "Must be dereferenceable to an integral type.");
        {
            std::string v;
            for(T i = it1; i < it2; v += std::to_string(*i++), v += ',');
            v.back() = '\n';
            ::std::cerr << "ks: " << v;
        }
        while(it1 < it2) {
            std::fprintf(stderr, "About to perform result batch for k = %d\n", *it1);
            ResultBatch rb(std::move(run(*it1++)));
            func(rb);
            std::fprintf(stderr, "Called function %s with k = %d\n", __PRETTY_FUNCTION__, *(it1 - 1));
        }
    }
};
FullAnalysis::ResultBatch::~ResultBatch() {}
FullAnalysis::~FullAnalysis() {}

void usage() {
    std::fprintf(stderr, "Usage: <executable> pair_paths.txt\nFlags:\n"
                         "-l: set linspace (required)\n"
                         "-k: add kmer size\n"
                         "-r: set kmer size range (int1,int2) (overrides previous -k)\n"
                         "-s: set sketch size range (int1,int2) [e.g., (16, 31)]\n"
                         "-m: use mash distance for bucketizing rather than jaccard index\n"
                         "-b: calculate bindash only\n"
                         "-o: set output file (/dev/stdout)\n"
                );
    std::exit(1);
}

ks::string make_fname(const char *fname, unsigned hllp, unsigned k, const std::string &scratch_folder, const char *suffix) {
    std::string pref = scratch_folder.size() ? std::string(scratch_folder + '/'): std::string("");
    std::string ps;
    const char *p = std::strrchr(fname, '/');
    if(p) ps = (p + 1);
    else ps = fname;
    return ks::sprintf("%s%s__.s%zu.k%u.%s", pref.data(), ps.data(), size_t(1) << (hllp - (2 + (k > 16))), k, suffix);
}

size_t bdnelem(size_t nbytes, unsigned bbits) {return (nbytes >> 3) / bbits / 4;}

static unsigned MAX_NTRIES = 20;

ks::string sketch_bindash(const char *path, unsigned hllp, unsigned k, const std::string &scratch_folder) {
    const size_t nels = bdnelem(size_t(1) << (hllp), 16); // Divide the usage by 8 (8 bytes in a u64), then 64 (64 u64s per nel), then 4 (I don't know why, but the sketch sizes were bigger than I expected.)
    auto outpath = make_fname(path, hllp, k, scratch_folder, "bdsh");
    if(isfile(outpath))
         return outpath;
    ks::string cstr = ks::sprintf("bindash sketch --minhashtype=2 --bbits=16 --sketchsize64=%zu --kmerlen=%u --outfname=%s %s", nels, k, outpath.data(), path);
    std::fprintf(stderr, "bindash command: %s\n", cstr.data());
    std::FILE *sketch;
    unsigned ntries = 0;
    do {
        if((sketch = popen(cstr.data(), "r")) == nullptr) {
            if(ntries++ > MAX_NTRIES) {
                std::fprintf(stderr, "[%s:%d:%s] Could not perform %s even after %u tries. errcode: %d/%s\n", __FILE__, __LINE__, __PRETTY_FUNCTION__, cstr.data(), ntries, errno, std::strerror(errno));
                std::fflush(stderr);
                std::exit(1);
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(1000));
            std::fprintf(stderr, "[%s:%s:%d] Failed to perform bindash command on try %u. Error code: %i. Message: %s\n", __FILE__, __PRETTY_FUNCTION__, __LINE__, ntries, errno, std::strerror(errno));
        }
    } while(sketch == nullptr);
    int retcode = pclose(sketch) >> 8;
    if(retcode) {
        std::fprintf(stderr, "sketch failed with error code %d for command %s. errno: %d/%s\n", retcode, cstr.data(), errno, std::strerror(errno));
        std::fflush(stderr);
        std::exit(1);
    }
    LOG_DEBUG("About to return mash sketch at %s\n", outpath.data());
    return outpath;
}

ks::string sketch_mash(const char *path, unsigned hllp, unsigned k, const std::string &scratch_folder) {
    const size_t nels = size_t(1) << (hllp >= (2 + 2 + (k > 16)) ? (hllp - (2 + 2 + (k > 16))): 0); // Divide the usage by 8 if k > 16, otherwise 4. (Also divide by 2 again, because the file sizes are bigger than expected.)
    auto outpath = make_fname(path, hllp, k, scratch_folder);
    if(isfile(outpath)) return outpath;
    ks::string cstr = ks::sprintf(MASHPATH " sketch -s %zu -k %u -o %s %s", nels, k, outpath.data(), path);
    std::FILE *sketch;
    unsigned ntries = 0;
    do {
        if((sketch = popen(cstr.data(), "r")) == nullptr) {
            if(ntries++ > MAX_NTRIES) {
                std::fprintf(stderr, "[%s:%d:%s] Could not perform %s even after %u tries. errcode: %d/%s\n", __FILE__, __LINE__, __PRETTY_FUNCTION__, cstr.data(), --ntries, errno, std::strerror(errno));
                std::fflush(stderr);
                std::exit(1);
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        }
        std::fprintf(stderr, "[%s:%s:%d] Failed to perform mash command on try %u. Error code: %i. Message: %s\n", __FILE__, __PRETTY_FUNCTION__, __LINE__, ntries, errno, std::strerror(errno));
    } while(sketch == nullptr);
    int retcode = pclose(sketch) >> 8;
    if(retcode) {
        std::fprintf(stderr, "sketch failed with error code %d for command %s. errno: %d/%s\n", retcode, cstr.data(), errno, std::strerror(errno));
        std::fflush(stderr);
        std::exit(1);
    }
    LOG_DEBUG("About to return mash sketch at %s\n", outpath.data());
    return outpath;
}

double get_bindash_ji(const char *path1, const char *path2) {
    std::FILE *fp;
    unsigned ntries = 0;
    const auto cmd = ks::sprintf("bindash dist --mthres=1e300 %s %s\n", path1, path2);
    LOG_INFO("Getting bindash ji with command %s\n", cmd.data());
    start:
    if((fp = popen(cmd.data(), "r")) == nullptr) {
        if(++ntries > MAX_NTRIES) throw std::runtime_error("Could not get bindash dist to run.");
        else goto start;
    }
    size_t ntabs = 0;
    int c;
    std::string s;
    for(int c;(c = std::fgetc(fp)) >= 0;) {
        s.push_back(c);
        if(c == '\t' && ++ntabs == 4) break;
    }
    if(ntabs != 4) {
        std::fprintf(stderr, "Sometimes bindash fails to report an answer. We return -1 for this case. We've found that re-performing the comparison does not enable its completion.");
        throw std::runtime_error("Number of tabs found was"s + std::to_string(ntabs) + ". This answer is invalid. Command: " + ks::sprintf("bindash dist %s %s\n", path1, path2).data());
    }
    std::string buf;
    while((c = std::fgetc(fp)) >= 0) {buf += c;}
    if(pclose(fp)) throw std::runtime_error("Error closing bindash dist call");
    double num = std::atof(buf.data());
    return num / std::atof(std::strchr(buf.data(), '/') + 1);
}

double get_mash_ji(const char *path1, const char *path2) {
    std::FILE *fp;
    std::string output;
    int ntries = 0;
    if(!isfile(path1) || !isfile(path2)) {
        LOG_INFO("Missing path or paths! %s, %s\n", path1, path2);
    }
    begin:
    if((fp = popen(ks::sprintf("mash dist -d 0. -j %s %s\n", path1, path2).data(), "r")) == nullptr) throw std::runtime_error("Could not get mash dist to run.");
    for(int c; (c = std::fgetc(fp)) >= 0; output += c);
    int rc;
    if((rc = pclose(fp)) != 0) {
        if(++ntries == 5)
            throw std::runtime_error("Error closing mash dist call");
        LOG_INFO("Mash dist call failed on try %i with retcode %d, errno %d and strerror %s. Try again...\n", ntries, rc, errno, strerror(errno));
        goto begin;
    }
    // I could actually parse this in a streaming fashion...
    char *p = &output[0];
    while(!std::isspace(*p)) ++p;
    while(std::isspace(*p)) ++p;
    while(!std::isspace(*p)) ++p;
    while(std::isspace(*p)) ++p;
    char *q = p + 1;
    while(!std::isspace(*q)) ++q;
    *q = '\0';
    return std::atof(p);
}

auto main(int argc, char *argv[]) -> int {
    std::FILE *fp = stdout;
    std::vector<int> kmersizes;
    int c;
    bool use_ji = true, bdonly = false;
    int pstart = 10, pend = 18, linspace = -1, nthreads = 1;
    std::string reference_ani_results, scratch_folder;
    while((c = getopt(argc, argv, "S:p:o:l:s:r:k:a:bh?m")) >= 0) {
        switch(c) {
            case 'b': bdonly = true; break;
            case 'p': nthreads = std::atoi(optarg); break;
            case 'a': reference_ani_results = optarg; break;
            case 'k': kmersizes.push_back(std::atoi(optarg)); break;
            case 'm': use_ji = false; break;
            case 'S': scratch_folder = optarg; break;
            case 'r': {
                    std::string arg = optarg;
                    if(arg.find(',') == arg.npos) {
                        kmersizes.push_back(std::atoi(arg.data()));
                        LOG_INFO("Note: only one kmer size provided.\n");
                        break;
                    }
                    char *p1 = &arg[0], *p2 = std::strchr(p1, ',');
                    *p2++ = '\0';
                    kmersizes.clear();
                    for(int i = std::atoi(p1); i <= std::atoi(p2); kmersizes.push_back(i++));
                    break;
            }
            case 'o': fp = std::fopen(optarg, "wb"); break;
            case 's': {
                    std::string arg = optarg;
                    if(arg.find(',') == arg.npos) {
                        pstart = pend = std::atoi(arg.data());
                        LOG_INFO("Note: only one sketch size provided.\n");
                        break;
                    }
                    char *p1 = &arg[0], *p2 = std::strchr(p1, ',');
                    *p2++ = '\0';
                    pstart = std::atoi(p1), pend = std::atoi(p2);
                    break;
            }
            case 'l': linspace = std::atoi(optarg); break;
        }
    }
    kmersizes = std::vector<int>{16, 21, 31};
    nthreads = std::max(nthreads, 1);
    omp_set_num_threads(nthreads);
    if(linspace < 1) usage();
    const char *inpath = argc == optind ? "": (const char *)argv[optind];
    if(pstart < 8 || pend < pstart || pend > 32) throw std::runtime_error("Nonsensical -s argument.");
    if(kmersizes.empty()) kmersizes = std::vector<int>{16, 21, 32};
    std::fprintf(fp, "#Path1\tPath2\tKmerSize\tBytes in Sketch\tExact\tOriginal\tMLE\tJMLE\tmash\n");
    std::fflush(fp);
    FullAnalysis analysis(0, inpath, pstart, pend, linspace, use_ji, scratch_folder, bdonly);
    analysis.krange_for_each(kmersizes.begin(), kmersizes.end(), [fp](const FullAnalysis::ResultBatch &rb){rb.write(fp); std::fflush(fp);});
}

