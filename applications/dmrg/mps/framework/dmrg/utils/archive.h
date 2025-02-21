/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2011-2012 by Michele Dolfi <dolfim@phys.ethz.ch>
 * 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#ifndef STORAGE_ARCHIVE_H
#define STORAGE_ARCHIVE_H

#include <boost/utility.hpp>
#include <alps/hdf5.hpp>
#include <alps/utility/encode.hpp>

namespace storage {

    inline std::string once(std::string fp){
        return fp;
    }

    inline void uniq(std::string fp){
    }

    class archive : boost::noncopyable {
    public:
        archive(std::string fp) : write(false), fp(fp) {
            impl = new alps::hdf5::archive(fp);
        }
        archive(std::string fp, const char* rights) : write(strcmp(rights,"w") == 0), fp(fp) {
            impl = new alps::hdf5::archive(once(fp), rights); 
        }
       ~archive(){
           delete impl;
           if(write) uniq(fp); 
        }
        bool is_group(const char* path){
            return impl->is_group(path);
        }
        bool is_scalar(const char* path){
            return impl->is_scalar(path);
        }
        bool is_data(const char* path){
            return impl->is_data(path);
        }
        template<typename T>
        void operator << (const T& obj){
            (*impl) << obj;
        }
        template<typename T>
        void operator >> (T& obj){
            (*impl) >> obj;
        }
        alps::hdf5::detail::archive_proxy<alps::hdf5::archive> operator[](std::string path){
            return (*impl)[path];
        }
    private:
        bool write;
        std::string fp;
        alps::hdf5::archive* impl;
    };
    
    inline std::string encode(std::string const & s){
        return alps::hdf5_name_encode(s);
    }
}

#endif
