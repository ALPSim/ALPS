/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005 by Lukas Gamper <mistral@student.ethz.ch>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
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
