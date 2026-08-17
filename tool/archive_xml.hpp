/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005 by Lukas Gamper <mistral@student.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
