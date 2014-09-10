#ifndef _MCPELE_SERIALIZE_VECTOR_H
#define _MCPELE_SERIALIZE_VECTOR_H

#include <stdio.h>
#include <sqlite3.h> //requires sudo apt-get install libsqlite3-dev
#include <iostream>
#include <sstream>
#include <cstdint>


namespace mcpele{

    /*Serialize and deserialize functions. These only work for vectors of Plain Old Data structures.
     * Furthermore these only work if the vectors contain only a single type and no pointers or references.
     * For compatibility reasons across platforms when storing arrays of integers must specify
     *  int32_t or uint32_t
    */

    //implementation for floating point values
    template<typename T>
    std::string serialize_vector(const std::vector<T>& vec)
    {
        static_assert(std::is_floating_point<T>::value && !std::is_compound<T>::value,"type is not floating point or is compound");
        std::ostringstream strm;
        strm.write(reinterpret_cast<const char*>(&vec[0]), vec.size()*sizeof(T));

        if (strm.fail()){
            strm.clear();
            throw std::runtime_error("binary file writing error");
        }

        return strm.str();
    }

    //implementation for fixed width integer 32 bit (for portability across platforms)
    template<>
    std::string serialize_vector(const std::vector<int32_t>& vec)
    {
        std::ostringstream strm;
        strm.write(reinterpret_cast<const char*>(&vec[0]), vec.size()*sizeof(int32_t));

        if (strm.fail()){
            strm.clear();
            throw std::runtime_error("binary file writing error");
        }

        return strm.str();
    }

    //implementation for fixed width unsigned integer 32 bit (for portability across platforms)
    template<>
    std::string serialize_vector(const std::vector<uint32_t>& vec)
    {
        std::ostringstream strm;
        strm.write(reinterpret_cast<const char*>(&vec[0]), vec.size()*sizeof(uint32_t));

        if (strm.fail()){
            strm.clear();
            throw std::runtime_error("binary file writing error");
        }

        return strm.str();
    }

    //implementation for floating point values
    template<typename T>
    std::vector<T> deserialize_vector(const std::string& ser_vec)
    {
        static_assert(std::is_floating_point<T>::value && !std::is_compound<T>::value,"type is not floating point or is compound");
        std::istringstream strm(ser_vec, std::iostream::binary);
        if (strm.fail()){
            strm.clear();
            throw std::runtime_error("binary file reading error");
        }
        const size_t length = ser_vec.size()/sizeof(T);
        std::vector<T> vec(length);
        strm.read(reinterpret_cast<char*>(&vec[0]), ser_vec.size());
        return vec;
    }

    //implementation for fixed width integer 32 bit (for portability across platforms)
    template<>
    std::vector<int32_t> deserialize_vector(const std::string& ser_vec)
    {
        std::istringstream strm(ser_vec, std::iostream::binary);
        if (strm.fail()){
            strm.clear();
            throw std::runtime_error("binary file reading error");
        }
        const size_t length = ser_vec.size()/sizeof(int32_t);
        std::vector<int32_t> vec(length);
        strm.read(reinterpret_cast<char*>(&vec[0]), ser_vec.size());
        return vec;
    }

    //implementation for fixed width unsigned integer 32 bit (for portability across platforms)
    template<>
    std::vector<uint32_t> deserialize_vector(const std::string& ser_vec)
    {
        std::istringstream strm(ser_vec, std::iostream::binary);
        if (strm.fail()){
            strm.clear();
            throw std::runtime_error("binary file reading error");
        }
        const size_t length = ser_vec.size()/sizeof(uint32_t);
        std::vector<uint32_t> vec(length);
        strm.read(reinterpret_cast<char*>(&vec[0]), ser_vec.size());
        return vec;
    }
}
#endif



