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


bool operator<(const std::vector<size_t>& lhs, const std::vector<size_t>& rhs) {
    for (size_t i = 0; i != lhs.size(); ++i) {
        if (lhs[i] < rhs[i]) {
            return true;
        }
    }
    return false;
}

class Bin {
public:
    Bin() : residuals_(nullptr), bitmap_(nullptr) {}
public:
    void unmapping(size_t ress, size_t sz) {
        munmap(residuals_, ress * sizeof(float));
	size_t tmp = sz / 8;
	if (sz % 8) {
            ++tmp;
	}
	munmap(bitmap_, tmp);
    }

    void* residuals_;
    void* bitmap_;
};


class Chunk {
public:
    Chunk() {}
    
    Chunk(const Chunk& other) : sz(other.sz), res_name(other.res_name), pos_name(other.pos_name) {}

    ~Chunk() {
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
        fdpos_ = open(chunk_.pos_name.data(), O_RDONLY);
        size_t tmp = sz / 8;
        if (sz % 8) {
            ++tmp;
        }
        tot_ops = tmp / sizeof(unsigned);
        if (tmp % sizeof(unsigned)) {
            ++tot_ops;
        }
        tmp = tot_ops * sizeof(unsigned);
        val_arr = mmap(0, chunk_.sz * sizeof(float), PROT_READ, MAP_SHARED, fdres_, 0);
	    LoadInfo();
    }

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
		    size_t tmp = chunk_sz_ / 8;
		    if (chunk_sz_ % 8) {
		        ++tmp;
		    }
            munmap(val_arr, chunk_.sz * sizeof(float));
            close(fdres_);
            close(fdpos_);
            return;
        }
        std::vector<size_t> offs;
        size_t ops = chunk_sz_ / (sizeof(unsigned) * 8);
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
        start_ = static_cast<float*>(mmap(0, full_sz, PROT_READ, MAP_SHARED, fd, 0));
        std::cerr << fd << " " << full_sz << " " << start_ << std::endl;
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
         std::cout << vec[i] << " ";
     }
     std::cout << "\n";
}

static int initial = 0;

class Database {
public:

    Database(std::map<std::vector<size_t>, Chunk>& tree) {
        tree_.swap(tree);
    }
    
    Database() = default;
    
    Database(Database&& other) {
        tree_.swap(other.tree_);
    }
    
    void simplend(const char* fname, std::vector<size_t> dims, std::vector<size_t> ch_sz, float lb, float ub, size_t bins, size_t bins_per_chunk) {
        ++initial;
        std::string init_str = std::to_string(initial);
        magic_number = bins_per_chunk;
        FileAccess file(fname, dims);
        std::vector<float> lbs;
        float dif = (ub - lb) / bins;
        for (size_t i = 0; i != bins; ++i) {
            lbs.push_back(lb + i * dif);
        }
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
		    if (act_sz % 8) {
		        act_sz /=  8;
		        ++act_sz;
		    } else {
		        act_sz /= 8;
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
				std::string namepos = maname + "pos" + init_str;
				int fdpos = open(namepos.data(), O_RDWR | O_CREAT | O_TRUNC, 0777);
				for (size_t j = 0; j != bins_per_chunk; ++j) {
				    cur_poss = std::vector<unsigned>(act_sz / sizeof(unsigned), 0);
				    write(fdres, ress[i * bins_per_chunk + j].data(), ress[i * bins_per_chunk + j].size() * sizeof(float));
				    total_vals += ress[i * bins_per_chunk + j].size();
				    for (size_t k = 0; k != poss[i * bins_per_chunk + j].size(); ++k) {
				        size_t cur = poss[i * bins_per_chunk + j][k];
				        size_t glob = cur / (sizeof(unsigned) * 8);
				        size_t loc = cur % (sizeof(unsigned) * 8);
				        unsigned sth = cur_poss[glob];
				        sth |= static_cast<unsigned>(1) << loc;
				        cur_poss[glob] = sth;
				    }
				    write(fdpos, cur_poss.data(), act_sz);
				}
				chunk.res_name = nameres;
				chunk.pos_name = namepos;
				chunk.sz = total_vals;
				tree_[coords] = chunk;
				chunk.res_name = "";
				chunk.pos_name = "";
				close(fdres);
				close(fdpos);
			}
        }  while (iter.next());
    }

    
public:
    std::map<std::vector<size_t>, Chunk> tree_;
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
    std::map<std::vector<size_t>, Chunk> res_tree;
    ++initial;
    for (auto it = lhs.tree_.begin(), jt = rhs.tree_.begin(); it != lhs.tree_.end(); ++it, ++jt) {
        auto& [lvec, lch] = *it;
        auto& [rvec, rch] = *jt;
        ChunkIterator binl(&lch, lvec[lvec.size() - 1]);
        ChunkIterator binr(&rch, rvec[rvec.size() - 1]);
        std::string maname = get_name(lvec);
        std::string init_str = std::to_string(initial);
		std::string nameres = maname + "res" + init_str;
		int fdres = open(nameres.data(), O_RDWR | O_CREAT | O_TRUNC, 0777);
		std::string namepos = maname + "pos" + init_str;
		int fdpos = open(namepos.data(), O_RDWR | O_CREAT | O_TRUNC, 0777);
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
                if (*(loffsets + i) == *(roffsets + j)) {
                    // std::cerr << lvalues << " " << rvalues << " " << lhs_sz << " " << i << " " << rhs_sz << " " << j << "\n";
                    if (lvalues[i] == rvalues[j]) {
                        poss.push_back(loffsets[i]);
                        ress.push_back(lvalues[i]);
                    }
                    ++i;
                    ++j;
                    continue;
                }
                if (*(loffsets + i) < *(roffsets + j)) {
                    ++i;
                } else {
                    ++j;
                }
            }
            ///////////////////////////////////////////////////////////
            size_t sz = lvec[lvec.size() - 1];
		    size_t act_sz = sz;
		    if (act_sz % 8) {
		        act_sz /=  8;
		        ++act_sz;
		    } else {
		        act_sz /= 8;
		    }
		    if (act_sz % sizeof(unsigned)) {
                act_sz = (act_sz / sizeof(unsigned) + 1) * sizeof(unsigned);
		    };
		    std::vector<unsigned> cur_poss = std::vector<unsigned>(act_sz / sizeof(unsigned), 0);
		    write(fdres, ress.data(), ress.size() * sizeof(float));
		    total_vals += ress.size();
		    for (size_t k = 0; k != poss.size(); ++k) {
		        size_t cur = poss[k];
		        size_t glob = cur / (sizeof(unsigned) * 8);
		        size_t loc = cur % (sizeof(unsigned) * 8);
		        unsigned sth = cur_poss[glob];
		        sth |= static_cast<unsigned>(1) << loc;
		        cur_poss[glob] = sth;
		    }
		    write(fdpos, cur_poss.data(), act_sz);           
            //////////////////////////////////////////////////////
        } while (binl.next_bin() && binr.next_bin());
		chunk.res_name = nameres;
		chunk.pos_name = namepos;
		chunk.sz = total_vals;
		close(fdres);
		close(fdpos);
        res_tree[lvec] = chunk;
        chunk.res_name = "";
        chunk.pos_name = "";
    }
    return Database(res_tree);
}




template <typename F>
Database value_similarity_join(Database& lhs, Database& rhs, F&& func, float tol) {
    std::map<std::vector<size_t>, Chunk> res_tree;
    for (auto it = lhs.tree_.begin(), jt = rhs.tree_.begin(); it != lhs.tree_.end(); ++it, ++jt) {
        auto& [lvec, lch] = *it;
        auto& [rvec, rch] = *jt;
        ChunkIterator binl(&lch, lvec[lvec.size() - 1]);
        ChunkIterator binr(&rch, rvec[rvec.size() - 1]);
        size_t bins_per_chunk = binl.bins_per_chunk_;
        Chunk chunk;
        Bin* bs = new Bin [bins_per_chunk];
        size_t cnt = 0;
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
            while ((i != lhs_sz) && (j != rhs_sz)) {
                if (*(loffsets + i) == *(roffsets + j)) {
                    if (func(lvalues[i], rvalues[j]) < tol) {
                        poss.push_back(loffsets[i]);
                        ress.push_back(lvalues[i]);
                    }
                    ++i;
                    continue;
                }
                if (*(loffsets + i) < *(roffsets + j)) {
                    ++i;
                } else {
                    ++j;
                }
            }
            ///////////////////////////////////////////////////////////
            bs[cnt].residuals_ = mmap(nullptr, ress.size() * sizeof(float), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, 0, 0);
            memcpy(bs[cnt].residuals_, ress.data(), ress.size() * sizeof(float));
			size_t sz = binl.chunk_sz_;
            if (sz % 8) {
                sz /=  8;
                ++sz;
            } else {
                sz /= 8;
            }
		    if (sz % sizeof(unsigned)) {
                sz = (sz / sizeof(unsigned) + 1) * sizeof(unsigned);
		    }
            unsigned* ptr = static_cast<unsigned*>(mmap(nullptr, sz, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, 0, 0));
            for (size_t k = 0; k * sizeof(unsigned) != sz; ++k) {
                *(ptr + k) = 0;
            }
            for (size_t k = 0; k != poss.size(); ++k) {
                size_t cur = poss[k];
                size_t glob = cur / (sizeof(unsigned) * 8);
                size_t loc = cur % (sizeof(unsigned) * 8);
                *(ptr + glob) |= (1 << loc);
            }
            bs[cnt].bitmap_ = ptr;
            ++cnt;
            //////////////////////////////////////////////////////
        } while (binl.next_bin() && binr.next_bin());
        chunk.bins_ = bs;
        res_tree[lvec] = chunk;
		chunk.bins_ = nullptr;
    }
    return Database(res_tree);
}











