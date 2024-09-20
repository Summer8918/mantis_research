#include <sdsl/vectors.hpp>
#include <sdsl/coder_elias_gamma.hpp>
#include <sdsl/coder_comma.hpp>
#include <sdsl/coder_fibonacci.hpp>
#include <iostream>
#include <sdsl/io.hpp>

using namespace sdsl;
using namespace std;

void test_int_vector() {
    std::string  tmp_file = "tmp_file.sdsl";
    {
        // generate int_vector and store to a file
        int_vector<> v = {1,3,5,7,2,3,4,9,8,7,10,1};
        store_to_file(v, tmp_file);
    }
    {
        // use int_vector_buffer to open the file
        int_vector_buffer<> ivb(tmp_file);
        // output elements and assign 42 to each of them
        for (auto x : ivb) {
            cout << x << endl;
            x = 42;
        }
        // int_vector_buffer is destroy at the end of the
        // scope and all value are written to disk
    }
    {
        // read vector from file and output it
        int_vector<> v;
        load_from_file(v, tmp_file);
        cout << v << endl;
    }
    // delete temporary file
    //sdsl::remove(tmp_file);
}

void test_int_vector_with_compression() {
    int_vector<32> v(10*(1<<20), 0);
    v[0] = 1ULL<<31;
    v[9 * (1<<20)] = 0xffffff;
    //util::bit_compress(v);
    std::string  tmp_file = "compressed_iv_tmp_file.sdsl";
    cout << "Uncompressed size in MB: " << size_in_mega_bytes(v) << endl;
    vlc_vector<coder::fibonacci> vv(v);
    cout << "Compressed size in MB: " << size_in_mega_bytes(vv) << endl;
	cout << "Percentage: " << size_in_mega_bytes(vv) / size_in_mega_bytes(v) * 100 << endl;
    store_to_file(vv, tmp_file);
    vlc_vector<coder::fibonacci> vv2;
    load_from_file(vv2, tmp_file);
    assert(vv2[0] == 1ULL<<63);
    assert(vv2[9 * (1<<20)] == 0xffffff);
    cout << vv2[0] << " " << vv2[9 * (1<<20)] << endl;
    cout << "int vector.capacity()" << v.capacity() << endl;
    cout << "int vector.max_size()" << v.max_size() << endl;
    cout << "int vector.bit_size()" << v.bit_size() << endl;
}

int main()
{
    test_int_vector();
    test_int_vector_with_compression();
    return 0;
}