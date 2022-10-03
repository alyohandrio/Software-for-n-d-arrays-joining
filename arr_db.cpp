#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <map>
#include <string.h>
#include <unistd.h>
#include <vector>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>


enum {
    BITS_PER_BYTE = 8,
};

struct VectorCmp {
    bool operator()(const std::vector<size_t>& lhs, const std::vector<size_t>& rhs) const {
        for (size_t i = 0; i != lhs.size(); ++i) {
            if (lhs[i] < rhs[i]) {
                return true;
            }
            if (lhs[i] > rhs[i]) {
                return false;
            }
        }
        return false;
    }
};

class Bin {
public:
    Bin() : residuals_(nullptr), bitmap_(nullptr) {}
public:
    void unmapping(size_t ress, size_t sz) {
        munmap(residuals_, ress * sizeof(float));
    size_t tmp = sz / BITS_PER_BYTE;
    if (sz % BITS_PER_BYTE) {
        ++tmp;
    }
        munmap(bitmap_, tmp);
    }

    void* residuals_;
    void* bitmap_;
};

// int alive = 0;

class Chunk {
public:
    Chunk() {/*++alive;*/}
    
    Chunk(const Chunk& other) : sz(other.sz), res_name(other.res_name), pos_name(other.pos_name) {/*++alive;*/}

    ~Chunk() {
        /*std::cerr << --alive << std::endl;*/
        if (res_name != "") {
            std::remove(res_name.data());
            std::remove(pos_name.data());
        }
    }

public:
    size_t sz;
    std::string res_name;
    std::string pos_name;
};


static size_t magic_number;


class ChunkIterator {
public:
    ChunkIterator(Chunk* chunk, size_t sz) : chunk_sz_(sz), chunk_(*chunk) {
        bins_per_chunk_ = magic_number;
        fdres_ = open(chunk_.res_name.data(), O_RDONLY);
        if (fdres_ < 0) {
            std::cerr << "ChunkIterator: Failed to open residuals file" << std::endl;
        }
        fdpos_ = open(chunk_.pos_name.data(), O_RDONLY);
        if (fdpos_ < 0) {
            std::cerr << "ChunkIterator: Failed to open positions file" << std::endl;
        }
        // std::cerr << "Opening " << fdres_ << ", " << fdpos_ << std::endl;
        size_t tmp = sz / BITS_PER_BYTE;
        if (sz % BITS_PER_BYTE) {
            ++tmp;
        }
        tot_ops = tmp / sizeof(unsigned);
        if (tmp % sizeof(unsigned)) {
            ++tot_ops;
        }
        // tmp = tot_ops * sizeof(unsigned);
        val_arr = mmap(0, chunk_.sz * sizeof(float), PROT_READ, MAP_SHARED, fdres_, 0);
        LoadInfo();
    }

    // ~ChunkIterator()

    void print() {
        do {
            std::cout << "Offsets:\n";
            for (size_t i = 0; i != offsets_.size(); ++i) {
                std::cout << offsets_[i] << " ";
            }
            std::cout << "\nValues:\n";
            for (size_t i = 0; i != offsets_.size(); ++i) {
                std::cout << *(values_ + i) << "  ";
            }
            std::cout << "\n\n";
        } while (next_bin());
    }
    
    bool next_bin() {
        if (pos == bins_per_chunk_) {
            return false;
        }
        ++pos;
        LoadInfo();
        return (pos != bins_per_chunk_);
    }

    size_t bin_size() {
        return offsets_.size();
    }
    
    size_t* offsets() {
        return offsets_.data();
    }
    
    float* values() {
        return values_;
    }
    
public:
    void LoadInfo() {
        if (pos == bins_per_chunk_) {
            /*size_t tmp = chunk_sz_ / BITS_PER_BYTE;
            if (chunk_sz_ % BITS_PER_BYTE) {
                ++tmp;
            }*/
            munmap(val_arr, chunk_.sz * sizeof(float));
            // std::cerr << "Closing " << fdres_ << ", " << fdpos_ << std::endl;
            close(fdres_);
            close(fdpos_);
            return;
        }
        std::vector<size_t> offs;
        size_t ops = chunk_sz_ / (sizeof(unsigned) * BITS_PER_BYTE);
        size_t cnt = 0;
        unsigned tmp;
        for (size_t i = 0; i != ops; ++i) {
            read(fdpos_, &tmp, sizeof(tmp));
            for (size_t j = 0; j != sizeof(unsigned); ++j) {
                if (tmp & (static_cast<unsigned>(1) << j)) {
                    offs.push_back(cnt);
                }
                ++cnt;
            }
        }
        ops = chunk_sz_ - ops * sizeof(unsigned) * 8;
        if (ops) {
            read(fdpos_, &tmp, sizeof(tmp));
        }
        for (size_t i = 0; i != ops; ++i) {
            if (tmp & (static_cast<unsigned>(1) << i)) {
                offs.push_back(cnt);
            }
            ++cnt;
        }
        offsets_.swap(offs);
        values_ = static_cast<float*>(val_arr) + floats_read;
        floats_read += offsets_.size();
    }


    size_t floats_read = 0;
    size_t tot_ops;
    size_t chunk_sz_;
    size_t bins_per_chunk_;
    Chunk& chunk_;
    size_t pos = 0;
    std::vector<size_t> offsets_;
    float* values_;
    int fdres_, fdpos_;
    void* val_arr;
};

class FileAccess {
public:
    FileAccess(const char* file, const std::vector<size_t>& dims) : dims_(dims) {
        full_sz = sizeof(float);
        size_t dims_sz = dims_.size();
        muls_.resize(dims_sz);
        for (size_t i = 0; i != dims_sz; ++i) {
            full_sz *= dims_[i];
        }
        muls_[dims_sz - 1] = 1;
        for (size_t i = 1; i != dims_sz; ++i) {
            muls_[dims_sz - 1 - i] = muls_[dims_sz - i] * dims_[dims_sz - i];
        }
        fd = open(file, O_RDONLY);
        if (fd < 0) {
            std::cerr << "Failed to open " << file << std::endl;
        }
        start_ = static_cast<float*>(mmap(NULL, full_sz, PROT_READ, MAP_SHARED, fd, 0));
        // std::cerr << "fd: " << fd << " full size: " << full_sz << " start: " << start_ << std::endl;
    }
    
    float get(std::vector<size_t>& coords) {
        size_t shift = 0;
        for (size_t i = 0; i != dims_.size(); ++i) {
            shift += coords[i] * muls_[i];
        }
        return *(start_ + shift);
    }
    
    ~FileAccess() {
        munmap(start_, full_sz);
        close(fd);
    }
private:
    std::vector<size_t> dims_;
    std::vector<size_t> muls_;
    float* start_;
    size_t full_sz;
    int fd;
};




class VectorIteration {
public:
    VectorIteration(const std::vector<size_t>& start, const std::vector<size_t>& fin) : it_(start), end_(fin), sz_(start.size()) {}
    
    bool next() {
        for (size_t i = 0; i != sz_; ++i) {
            if (++it_[sz_ - 1 - i] < end_[sz_ - i - 1]) {
                return true;
            }
            it_[sz_ - i - 1] = 0;
        }
        return false;
    }
    
    size_t operator[](size_t i) {
        return it_[i];
    }
    
    size_t pos() {
        size_t res = it_[0];
        for (size_t i = 1; i != sz_; ++i) {
            res *= end_[i];
            res += it_[i];
        }
        return res;
    }
    
    std::vector<size_t> vec() {
        return it_;
    }
    
private:
    std::vector<size_t> it_;
    std::vector<size_t> end_;
    size_t sz_;
};


void safe_write(int fd, char* data, size_t sz) {
    ssize_t wrote;
    size_t written = 0;
    while (written != sz) {
        wrote = write(fd, data + written, sz - written);
        if (wrote < 0) {
            std::cerr << "Error during safe_write occured with fd " << fd << std::endl;
            _exit(1);
        } else {
            written += wrote;
        }
    }
}


std::string get_name(const std::vector<size_t>& vec) {
    std::string res = "store";
    for (size_t i = 0; i != vec.size(); ++i) {
        res = res + "_" + std::to_string(vec[i]);
    }
    return res;
}


template <typename T>
void print_vec(const std::vector<T>& vec) {
     for (size_t i = 0; i != vec.size(); ++i) {
         std::cerr << vec[i] << " ";
     }
     std::cerr << std::endl;
}

static int initial = 0;

class Database {
public:

    Database(std::map<std::vector<size_t>, Chunk, VectorCmp>& tree, std::vector<size_t>& dims, std::vector<size_t>& ch_sz, std::vector<float>& lbs, float step, size_t bpc) :
     dims_(dims), ch_sz_(ch_sz), lbs_(lbs), step_(step), bins_per_chunk_(bpc) {
        tree_.swap(tree);
    }
    
    Database() = default;
    
    Database(Database&& other) = default;
    
    void simplend(const char* fname, std::vector<size_t> dims, std::vector<size_t> ch_sz, float lb, float ub, size_t bins, size_t bins_per_chunk) {
        dims_ = dims;
        ch_sz_ = ch_sz;
        bins_per_chunk_ = bins_per_chunk;
        ++initial;
        std::string init_str = std::to_string(initial);
        magic_number = bins_per_chunk;
        FileAccess file(fname, dims);
        std::vector<float> lbs;
        float dif = (ub - lb) / bins;
        for (size_t i = 0; i != bins; ++i) {
            lbs.push_back(lb + i * dif);
        }
        lbs_ = lbs;
        step_ = dif;
        std::vector<size_t> chunks_num(dims.size());
        for (size_t i = 0; i != dims.size(); ++i) {
            chunks_num[i] = dims[i] / ch_sz[i];
            if (dims[i] % ch_sz[i]) {
                ++chunks_num[i];
            }
        }
        std::vector<size_t> zeros(dims.size(), 0);
        VectorIteration iter(zeros, chunks_num);
        do {
            std::vector<std::vector<float>> ress(bins);
            std::vector<std::vector<size_t>> poss(bins);
            std::vector<size_t> max_bound = ch_sz;
            size_t sz = 1;
            for (size_t i = 0; i != dims.size(); ++i) {
                if ((iter[i] + 1) * ch_sz[i] > dims[i]) {
                    max_bound[i] = dims[i] - iter[i] * ch_sz[i];
                }
                sz *= max_bound[i];
            }
            VectorIteration jter(zeros, max_bound);
            do {
                size_t k = 0;
                std::vector<size_t> global_coords(dims.size());
                for (size_t i = 0; i != dims.size(); ++i) {
                    global_coords[i] = iter[i] * ch_sz[i] + jter[i];
                }
                float cur = file.get(global_coords);
                while ((k + 1 < bins) && (cur >= lbs[k + 1])) {
                    ++k;
                }
                ress[k].push_back(cur - lbs[k]);
                poss[k].push_back(jter.pos());
            } while (jter.next());
            size_t act_sz = sz;
            if (act_sz % BITS_PER_BYTE) {
                act_sz /= BITS_PER_BYTE;
                ++act_sz;
            } else {
                act_sz /= BITS_PER_BYTE;
            }
            if (act_sz % sizeof(unsigned)) {
                act_sz = (act_sz / sizeof(unsigned) + 1) * sizeof(unsigned);
            }
            std::vector<unsigned> cur_poss(act_sz / sizeof(unsigned), 0);
            for (size_t i = 0; i * bins_per_chunk < bins; ++i) {
                Chunk chunk;
                size_t total_vals = 0;
                std::vector<size_t> coords = iter.vec();
                coords.push_back(i);
                coords.push_back(sz);
                std::string maname = get_name(coords);
                std::string nameres = maname + "res" + init_str;
                int fdres = open(nameres.data(), O_RDWR | O_CREAT | O_TRUNC, 0777);
                if (fdres < 0) {
                    std::cerr << "Simplend: Failed to open residuals file" << std::endl;
                }
                std::string namepos = maname + "pos" + init_str;
                int fdpos = open(namepos.data(), O_RDWR | O_CREAT | O_TRUNC, 0777);
                if (fdres < 0) {
                    std::cerr << "Simplend: Failed to open positions file" << std::endl;
                }
                for (size_t j = 0; j != bins_per_chunk; ++j) {
                    cur_poss = std::vector<unsigned>(act_sz / sizeof(unsigned), 0);
                    safe_write(fdres, reinterpret_cast<char*>(ress[i * bins_per_chunk + j].data()), ress[i * bins_per_chunk + j].size() * sizeof(float));
                    total_vals += ress[i * bins_per_chunk + j].size();
                    for (size_t k = 0; k != poss[i * bins_per_chunk + j].size(); ++k) {
                        size_t cur = poss[i * bins_per_chunk + j][k];
                        size_t glob = cur / (sizeof(unsigned) * BITS_PER_BYTE);
                        size_t loc = cur % (sizeof(unsigned) * BITS_PER_BYTE);
                        unsigned sth = cur_poss[glob];
                        sth |= static_cast<unsigned>(1) << loc;
                        cur_poss[glob] = sth;
                    }
                    safe_write(fdpos, reinterpret_cast<char*>(cur_poss.data()), act_sz);
                }
                chunk.res_name = nameres;
                chunk.pos_name = namepos;
                chunk.sz = total_vals;
                tree_[coords] = chunk;
                chunk.res_name = ""; // This object will be destroyed now
                chunk.pos_name = "";
                close(fdres);
                close(fdpos);
            }
        }  while (iter.next());
    }


    /*void write_to_file(std::string filename) {
        size_t sz = 1;
        for (int dim : dims_) {
            sz *= dim;
        }
        int fd = open(filename.data(), O_RDWR | O_CREAT | O_TRUNC, 0777);
        if (ftruncate(fd, sz * sizeof(float)) < 0) {
            std::cerr << "Failed to write to file" << std::endl;
            return;
        }
        float *arr = mmap(NULL, sz * sizeof(float), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        if (!arr) {
            std::cerr << "Failed to write to file" << std::endl;
            return;
        }
        for (size_t i = 0; i != sz; ++i) {
            arr[i] = std::numeric_limits<float>::quiet_NaN();
        }
        
    }*/

    
public:
    std::map<std::vector<size_t>, Chunk, VectorCmp> tree_;
    std::vector<size_t> dims_;
    std::vector<size_t> ch_sz_;
    std::vector<float> lbs_;
    float step_;
    size_t bins_per_chunk_;
};



class DatabaseIterator {
public:
    DatabaseIterator(Database* db) : db_(*db) {
        it_ = db_.tree_.begin();
    }


    bool next() {
        if (it_ == db_.tree_.end()) {
            return false;
        }
        ++it_;
        return (it_ != db_.tree_.end());
    }

    void print() {
        do {
            auto& [vec, ch] = *it_;
            /*std::cout << "Chunk ( ";
            for (size_t i = 0; i + 2 != vec.size(); ++i) {
                std::cout << vec[i] << " ";
            }
            std::cout << "), bins " << vec[vec.size() - 2] << std::endl;*/
            ChunkIterator chit(&ch, vec[vec.size() - 1]);
            chit.print();
            puts("---------------");
        } while (next());
    }
    
    bool seek(std::vector<size_t>& id) {
        it_ = db_.tree_.find(id);
        return (it_ != db_.tree_.end());
    }
public:
    Database& db_;
    std::map<std::vector<size_t>, Chunk>::iterator it_;
};


Database equal_join(Database& lhs, Database& rhs) {
    std::map<std::vector<size_t>, Chunk, VectorCmp> res_tree;
    ++initial;
    for (auto it = lhs.tree_.begin(), jt = rhs.tree_.begin(); it != lhs.tree_.end(); ++it, ++jt) {
        auto& [lvec, lchunk] = *it;
        auto& [rvec, rchunk] = *jt;
        ChunkIterator binl(&lchunk, lvec[lvec.size() - 1]);
        ChunkIterator binr(&rchunk, rvec[rvec.size() - 1]);
        std::string maname = get_name(lvec);
        std::string init_str = std::to_string(initial);
        std::string nameres = maname + "res" + init_str;
        int fdres = open(nameres.data(), O_RDWR | O_CREAT | O_TRUNC, 0777);
        // std::cerr << "Res " << fdres << std::endl;
        if (fdres < 0) {
            std::cerr << "EqualJoin: Failed to open residuals file" << std::endl;
        }
        std::string namepos = maname + "pos" + init_str;
        int fdpos = open(namepos.data(), O_RDWR | O_CREAT | O_TRUNC, 0777);
        // std::cerr << "Pos " << fdpos << std::endl;
        if (fdpos < 0) {
            std::cerr << "EqualJoin: Failed to open positions file" << std::endl;
        }
        Chunk chunk;
        size_t total_vals = 0;
        do {
            std::vector<size_t> poss;
            std::vector<float> ress;
            size_t lhs_sz = binl.bin_size();
            size_t rhs_sz = binr.bin_size();
            float* lvalues = binl.values();
            float* rvalues = binr.values();
            size_t* loffsets = binl.offsets();
            size_t* roffsets = binr.offsets();
            size_t i = 0, j = 0;
            while ((i < lhs_sz) && (j < rhs_sz)) {
                if (loffsets[i] == roffsets[j]) {
                    if (lvalues[i] == rvalues[j]) {
                        poss.push_back(loffsets[i]);
                        ress.push_back(lvalues[i]);
                    }
                    ++i;
                    ++j;
                    continue;
                }
                if (loffsets[i] < roffsets[j]) {
                    ++i;
                } else {
                    ++j;
                }
            }
            ///////////////////////////////////////////////////////////
            size_t sz = lvec[lvec.size() - 1];
            size_t act_sz = sz;
            if (act_sz % BITS_PER_BYTE) {
                act_sz /= BITS_PER_BYTE;
                ++act_sz;
            } else {
                act_sz /= BITS_PER_BYTE;
            }
            if (act_sz % sizeof(unsigned)) {
                act_sz = (act_sz / sizeof(unsigned) + 1) * sizeof(unsigned);
            };
            std::vector<unsigned> cur_poss = std::vector<unsigned>(act_sz / sizeof(unsigned), 0);
            safe_write(fdres, reinterpret_cast<char*>(ress.data()), ress.size() * sizeof(float));
            total_vals += ress.size();
            for (size_t k = 0; k != poss.size(); ++k) {
                size_t cur = poss[k];
                size_t glob = cur / (sizeof(unsigned) * BITS_PER_BYTE);
                size_t loc = cur % (sizeof(unsigned) * BITS_PER_BYTE);
                unsigned sth = cur_poss[glob];
                sth |= static_cast<unsigned>(1) << loc;
                cur_poss[glob] = sth;
            }
            safe_write(fdpos, reinterpret_cast<char*>(cur_poss.data()), act_sz);
            //////////////////////////////////////////////////////
        } while (binl.next_bin() && binr.next_bin());
        binr.next_bin();
        chunk.res_name = nameres;
        chunk.pos_name = namepos;
        chunk.sz = total_vals;
        close(fdres);
        close(fdpos);
        res_tree[lvec] = chunk;
        chunk.res_name = "";
        chunk.pos_name = "";
    }
    return Database(res_tree, lhs.dims_, lhs.ch_sz_, lhs.lbs_, lhs.step_, lhs.bins_per_chunk_);
}







template <typename F>
Database value_similarity_join(Database& lhs, Database& rhs, F&& func, float tol) {
    std::map<std::vector<size_t>, Chunk, VectorCmp> res_tree;
    ++initial;
    size_t current_bin;
    // size_t epoch = 0;
    for (auto it = lhs.tree_.begin(); it != lhs.tree_.end(); ++it) {
        auto& [lvec, lchunk] = *it;
        if (lvec[lvec.size() - 2] == 0) {
            current_bin = 0;
        }
        ChunkIterator binl(&lchunk, lvec[lvec.size() - 1]);
        std::string maname = get_name(lvec);
        std::string init_str = std::to_string(initial);
        std::string nameres = maname + "res" + init_str;
        int fdres = open(nameres.data(), O_RDWR | O_CREAT | O_TRUNC, 0777);
        if (fdres < 0) {
            std::cerr << "EqualJoin: Failed to open residuals file" << std::endl;
        }
        std::string namepos = maname + "pos" + init_str;
        int fdpos = open(namepos.data(), O_RDWR | O_CREAT | O_TRUNC, 0777);
        if (fdpos < 0) {
            std::cerr << "EqualJoin: Failed to open positions file" << std::endl;
        }
        Chunk chunk;
        size_t total_vals = 0;
        std::vector<size_t> coords = lvec;
        coords[coords.size() - 2] = 0;
        auto jt = rhs.tree_.find(coords);
        size_t current_bin_r = 0;
        size_t start_current_bin = current_bin;
        // std::cerr << epoch++ << std::endl;
        while (jt != rhs.tree_.end()) {
            current_bin = start_current_bin;
            // print_vec(coords);
            // std::cerr << current_bin << '\t' << current_bin_r << std::endl;
            if ((func(lhs.lbs_[current_bin], rhs.lbs_[std::min(current_bin_r + rhs.bins_per_chunk_, rhs.lbs_.size()) - 1] + rhs.step_) >= tol) &&\
                 (func(lhs.lbs_[std::min(current_bin + lhs.bins_per_chunk_, lhs.lbs_.size()) - 1] + lhs.step_, rhs.lbs_[current_bin_r]) >= tol)) {
                current_bin_r += rhs.bins_per_chunk_;
                ++coords[coords.size() - 2];
                jt = rhs.tree_.find(coords);
                continue;
            }
            auto& [rvec, rchunk] = *jt;
            size_t start_current_bin_r = current_bin_r;
            /*if (start_current_bin_r > 0) {
                std::cerr << start_current_bin_r << std::endl;
            }*/
            do {
                ChunkIterator binr(&rchunk, rvec[rvec.size() - 1]);
                struct PR {
                    PR(size_t x, float y) : pos(x), res(y) {}
                    
                    size_t pos;
                    float res;
                };
                std::vector<PR> prs;
                current_bin_r = start_current_bin_r;
                do {
                    // std::cerr << current_bin << '\t' << current_bin_r << std::endl;
                    if ((func(lhs.lbs_[current_bin], rhs.lbs_[current_bin_r] + rhs.step_) >= tol) && (func(lhs.lbs_[current_bin] + lhs.step_, rhs.lbs_[current_bin_r]) >= tol)) {
                        ++current_bin_r;
                        continue;
                    }
                    size_t lhs_sz = binl.bin_size();
                    size_t rhs_sz = binr.bin_size();
                    float* lvalues = binl.values();
                    float* rvalues = binr.values();
                    size_t* loffsets = binl.offsets();
                    size_t* roffsets = binr.offsets();
                    size_t i = 0, j = 0;
                    while ((i < lhs_sz) && (j < rhs_sz)) {
                        if (loffsets[i] == roffsets[j]) {
                            // std::cerr << lvalues << " " << rvalues << " " << lhs_sz << " " << i << " " << rhs_sz << " " << j << "\n";
                            // std::cerr << current_bin << '\t' << current_bin_r << std::endl;
                            if (std::forward<F>(func)(lhs.lbs_[current_bin] + lvalues[i], rhs.lbs_[current_bin_r] + rvalues[j]) < tol) {
                                prs.emplace_back(loffsets[i], lvalues[i]);
                            }
                            ++i;
                            ++j;
                            continue;
                        }
                        if (loffsets[i] < roffsets[j]) {
                            ++i;
                        } else {
                            ++j;
                        }
                    }
                    ///////////////////////////////////////////////////////////
                    // std::cerr << '-' << std::endl;
                    ++current_bin_r;
                    // std::cerr << "Bin " << current_bin << std::endl;
                    //////////////////////////////////////////////////////
                } while (binr.next_bin());
                size_t sz = lvec[lvec.size() - 1];
                size_t act_sz = sz;
                if (act_sz % BITS_PER_BYTE) {
                    act_sz /= BITS_PER_BYTE;
                    ++act_sz;
                } else {
                    act_sz /= BITS_PER_BYTE;
                }
                if (act_sz % sizeof(unsigned)) {
                    act_sz = (act_sz / sizeof(unsigned) + 1) * sizeof(unsigned);
                };
                std::vector<unsigned> cur_poss = std::vector<unsigned>(act_sz / sizeof(unsigned), 0);
                std::sort(prs.begin(), prs.end(), [](const PR& lhs, const PR& rhs) {return lhs.pos < rhs.pos;});
                std::vector<float> ress(prs.size());
                std::vector<size_t> poss(prs.size());
                for (size_t cnt = 0; cnt != prs.size(); ++cnt) {
                    ress[cnt] = prs[cnt].res;
                    poss[cnt] = prs[cnt].pos;
                }
                // std::cerr << "Val-sim: writing residuals with fd = " << fdres << std::endl;
                safe_write(fdres, reinterpret_cast<char*>(ress.data()), ress.size() * sizeof(float));
                // std::cerr << '-' << std::endl;
                total_vals += ress.size();
                for (size_t k = 0; k != poss.size(); ++k) {
                    size_t cur = poss[k];
                    size_t glob = cur / (sizeof(unsigned) * BITS_PER_BYTE);
                    size_t loc = cur % (sizeof(unsigned) * BITS_PER_BYTE);
                    unsigned sth = cur_poss[glob];
                    sth |= static_cast<unsigned>(1) << loc;
                    cur_poss[glob] = sth;
                }
                // std::cerr << "Val-sim: writing positions with fd = " << fdpos << std::endl;
                safe_write(fdpos, reinterpret_cast<char*>(cur_poss.data()), cur_poss.size() * sizeof(unsigned));
                ++current_bin;
            } while (binl.next_bin());
            ++coords[coords.size() - 2];
            jt = rhs.tree_.find(coords);
        }

        chunk.res_name = nameres;
        chunk.pos_name = namepos;
        chunk.sz = total_vals;
        close(fdres);
        close(fdpos);
        res_tree[lvec] = chunk;
        chunk.res_name = "";
        chunk.pos_name = "";
    }
    return Database(res_tree, lhs.dims_, lhs.ch_sz_, lhs.lbs_, lhs.step_, lhs.bins_per_chunk_);
}

