/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005 by Lukas Gamper <mistral@student.ethz.ch>
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

/* $Id$ */

#ifndef _XML_H_
#define _XML_H_

#include <string>

#include <boost/filesystem/path.hpp>

namespace fs = boost::filesystem;

#include "archive_node.hpp"

/**
 * Function-object to parse xml-documents in an object-tree. The document is 
 * represented by a tree of node-objects. Each element is represented by a node.
 * 
 * @see Node
 */
class XML {
        bool mVerbose;
    std::string readFile(fs::path filename);
        
        public:        
                /**
                 * Constructor of the Class
                 * 
                 * @param verbose Decides if the Actions should be printed out to std::cout
                 */
                XML(bool verbose): mVerbose(verbose) {}
                
                /**
                 * Functionoperator to parse a file, specified by the path to a object tree
                 * 
                 * @param inFileName Path to the File, that should be parsed
                 */
                Node operator()(fs::path inFileName, bool usePlotDTD);
                
};
#endif //_XML_H_
