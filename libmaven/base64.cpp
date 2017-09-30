#include "base64.h"
#include "mzUtils.h"

using namespace std;

namespace base64 { 

std::string decode_base64(const std::string& in) {
	std::string out;

	std::vector<int> T(256,-1);
	for (int i=0; i<64; i++) T["ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"[i]] = i; 

	int val=0, valb=-8;
	for (unsigned char c : in) {
		if (T[c] == -1) break;
		val = (val<<6) + T[c];
		valb += 6;
		if (valb>=0) {
			out.push_back(char((val>>valb)&0xFF));
			valb-=8;
		}
	}
	return out;
}


vector<float> decode_base64(const string& src, int float_size, bool neworkorder, bool decompress) {

	string destStr=decode_base64(src);

	if (decompress) {
#ifdef ZLIB
		destStr =mzUtils::decompress_string(destStr);
#endif 
	}



#if (LITTLE_ENDIAN == 1)
	 cerr << "WARNING: LITTLE_ENDIAN.. Inverted network order";
     neworkorder=!neworkorder;
#endif

    int size = 1+(destStr.size() * 3/4 - 4)/float_size;

    //we will cast everything as a float may be this is not wise, but have not found a need for double
    //precission yet
    vector<float> decodedArray(size);
	const char* dest = destStr.c_str();

    ////cerr << "Net=" << neworkorder << " float_size=" << float_size << endl;
    if ( float_size == 8 ) {
        if ( neworkorder == false ) {
            for (int i=0; i<size; i++) decodedArray[i] = (float) ((double*)dest)[i];
        } else {
            uint64_t *u = (uint64_t *) dest;
            double data=0;
            for (int i=0; i<size; i++) {
                uint64_t t = swapbytes64(u[i]);
                memcpy(&data,&t,8);
                decodedArray[i] = (float) data;
            }
        }
    } else if (float_size == 4 ) {
        if ( neworkorder == false) {
            for (int i=0; i<size; i++) decodedArray[i] = ((float*)dest)[i];
        } else {
            uint32_t *u = (uint32_t *) dest;
            float data=0;
            for (int i=0; i<size; i++) {
                uint32_t t = swapbytes(u[i]);
                memcpy(&data,&t,4);
                decodedArray[i] = data;
            }
        }
    } 
    //for debuging
    //for(int i=0; i < size; i++ ) { cout << "\tok.." << i << " " << size <<" " << decodedArray[i] << endl; }

    return decodedArray;

}

}//namespace
