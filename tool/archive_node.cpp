/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005-2008 by Lukas Gamper <mistral@student.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
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

#include <stdexcept>

#include "archive_node.hpp"

void Node::addElement(Node inChild) {
        inChild.mParent = this;
        mElem.push_back(inChild);
}

void Node::addElement(std::string inName) {
        std::transform (inName.begin(), inName.end(), inName.begin(), tolower);
        addElement(Node(inName, this));
}

void Node::addAttribute(std::string inName, std::string inValue) {
        std::transform (inName.begin(), inName.end(), inName.begin(), tolower);
        mAttr[inName] = inValue;
}

std::string Node::getAttribute(std::string inName) const {
  std::transform (inName.begin(), inName.end(), inName.begin(), tolower);
  std::map<std::string, std::string>::const_iterator itr = mAttr.find(inName);
  if (itr == mAttr.end())
    return "";
  else
    return itr->second;
}

std::string Node::string() {
        std::string value;
        for (std::list<std::string>::iterator it = mText.begin(); it != mText.end(); it++)
                value += *it;
        return value;
}

std::list<Node> Node::nodeTest(std::string nodeName) {
        std::transform (nodeName.begin(), nodeName.end(), nodeName.begin(), tolower);
        std::list<Node> context;
        for (std::list<Node>::iterator it = mElem.begin(); it != mElem.end(); it++)
                if (it->getName() == nodeName)
                        context.push_back(*it);
        return context;
}

Node Node::getElement(std::string nodeName) {
        std::transform (nodeName.begin(), nodeName.end(), nodeName.begin(), tolower);
        for (std::list<Node>::iterator it = mElem.begin(); it != mElem.end(); it++)
                if (it->getName() == nodeName)
                        return *it;
        return Node();
}
